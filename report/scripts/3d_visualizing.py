import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from scipy.interpolate import interp2d, griddata

def p(x, y):
    """Returns the [real, imaginary] part of the pressure function P(x,y,0,omega)
    """
    ans = np.zeros((2,) + x.shape)
    for x_S, y_S, z_S in sources:
        dr = np.sqrt((x-x_S)**2 + (y-y_S)**2 + z_S**2)
        ans[0] += np.cos(k * dr) / dr
        ans[1] += -np.sin(k * dr) / dr
    return ans

def dG(x, y, x_A, y_A):
    """Returns the [real, imaginary] part of the Green's function (dG/dz)_{z=0}
    """
    dr = np.sqrt((x-x_A)**2 + (y-y_A)**2 + z_A**2)
    real_part1 =  z_A  /(dr**3) * np.cos(k*dr)
    real_part2 =  z_A*k/(dr**2) * np.sin(k*dr)
    imag_part1 = -z_A  /(dr**3) * np.sin(k*dr)
    imag_part2 =  z_A*k/(dr**2) * np.cos(k*dr)
    return [real_part1+real_part2, imag_part1+imag_part2]

def prod(x, y, x_A, y_A):
    """Returns the [real, imaginary] part of the product of p and dG
    """
    pressure = p(x, y)
    greens = dG(x, y, x_A, y_A)
    return [pressure[0]*greens[0] - pressure[1]*greens[1],
            pressure[0]*greens[1] + pressure[1]*greens[0]]


def integral(x, y, res, method='trapezeoid'):
    """Calculate the integral using various methods
    """
    if method == 'trapezeoid':  # equidistant is assumed
        h_x = inter_x/(x.shape[1]-1)  # average
        h_y = inter_y/(x.shape[0]-1)  # average
        c_x = np.array([1] + [2]*(x.shape[0]-2) + [1]) * h_x/2
        c_y = np.array([1] + [2]*(x.shape[1]-2) + [1]) * h_y/2
        C = c_x.reshape((len(x), 1)) * c_y
        real = np.multiply(C, res[0]).sum()
        imag = np.multiply(C, res[1]).sum()
        return [real, imag]

    elif method == 'altered_trapezeoid':  # rectilinear is assumed
        C = np.zeros(x.shape)
        for i_x in range(x.shape[1]-1):
            for i_y in range(x.shape[0]-1):
                h_x = (x[i_y][i_x+1]-x[i_y][i_x])
                h_y = (y[i_y+1][i_x]-y[i_y][i_x])
                h = h_x * h_y / 4
                C[i_y:i_y+2, i_x:i_x+2] += h

        real = np.multiply(C, res[0]).sum()
        imag = np.multiply(C, res[1]).sum()
        return [real, imag]

    elif method == 'simpson':  # equidistant is assumed
        lst = [1, 4, 1]
        frac = 1/3

        ### check if dimensions add up, this would not be necessary
        #   if using appropriate boundary interpolations
        if ((x.shape[0]-len(lst)) % (len(lst)-1) != 0
            or (x.shape[1]-len(lst)) % (len(lst)-1) != 0):
            print('change #ticks to', (x.shape[0]-len(lst)) % (len(lst)-1)
                                      (x.shape[1]-len(lst)) % (len(lst)-1))
            raise BaseException

        ### do the calculation
        h_x = inter_x/(samples_x-1)  # average
        h_y = inter_y/(samples_y-1)  # average

        c_x = np.array(
            [lst[0]] +
            (lst[1:-1] + [lst[-1]+lst[0]]) * ((x.shape[0]-len(lst))//(len(lst)-1)) +
            lst[1:]
        ) * (h_x*frac)

        c_y = np.array(
            [lst[0]] +
            (lst[1:-1] + [lst[-1]+lst[0]]) * ((x.shape[1]-len(lst))//(len(lst)-1)) +
            lst[1:]
        ) * (h_y*frac)

        C = c_x.reshape((x.shape[0], 1)) * c_y
        real = np.multiply(C, res[0]).sum()
        imag = np.multiply(C, res[1]).sum()
        return [real, imag]

    elif method == 'altered_simpson':  # semi-equidistant is assumed
        ### do the calculation
        tot_real = 0
        tot_imag = 0
        # looping through all possible starting points
        for i_x in range(0, len(x[0])-2, 2):
            for i_y in range(0, len(x)-2, 2):
                for is_complex in range(0, 2):  # real, complex
                    f = {(i%3, i//3):
                            res[is_complex][i_y+i%3, i_x+i//3] for i in range(9)}

                    h_x = sum((x[i][i_x+2] - x[i][i_x])/2
                                for i in range(i_y, i_y+3)) / 3  # average
                    h_y = sum((y[i_y+2][i] - y[i_y][i])/2
                                for i in range(i_x, i_x+3)) / 3  # average

                    d_x0 = x[i_y+0][i_x+1] - (x[i_y+0][i_x+2]+x[i_y+0][i_x])/2
                    d_x1 = x[i_y+1][i_x+1] - (x[i_y+1][i_x+2]+x[i_y+1][i_x])/2
                    d_x2 = x[i_y+2][i_x+1] - (x[i_y+2][i_x+2]+x[i_y+2][i_x])/2

                    d_y0 = y[i_y+1][i_x+0] - (y[i_y+2][i_x+0]+y[i_y][i_x+0])/2
                    d_y1 = y[i_y+1][i_x+1] - (y[i_y+2][i_x+1]+y[i_y][i_x+1])/2
                    d_y2 = y[i_y+1][i_x+2] - (y[i_y+2][i_x+2]+y[i_y][i_x+2])/2
                    d_y = (d_y0+d_y1+d_y2)/3  # no rotation of axis for better result


                    a_0 = h_x * (
                            -3*d_x0**2*(f[(0,0)]+f[(2,0)]) +
                             2*d_x0   *(f[(0,0)]-f[(2,0)])*h_x +
                               h_x**2 *(f[(0,0)]+4*f[(1,0)]+f[(2,0)])
                    ) / (3*(h_x**2-d_x0**2))

                    a_1 = h_x * (
                            -3*d_x1**2*(f[(0,1)]+f[(2,1)]) +
                             2*d_x1   *(f[(0,1)]-f[(2,1)])*h_x +
                               h_x**2 *(f[(0,1)]+4*f[(1,1)]+f[(2,1)])
                    ) / (3*(h_x**2-d_x1**2))

                    a_2 = h_x * (
                            -3*d_x2**2*(f[(0,2)]+f[(2,2)]) +
                             2*d_x2   *(f[(0,2)]-f[(2,2)])*h_x +
                               h_x**2 *(f[(0,2)]+4*f[(1,2)]+f[(2,2)])
                    ) / (3*(h_x**2-d_x2**2))


                    if is_complex == 0:
                        tot_real += h_y * (
                                -3*d_y**2*(a_0+a_2) +
                                 2*d_y   *(a_0-a_2)*h_y +
                                   h_y**2*(a_0+4*a_1+a_2)
                        ) / (3*(h_y**2-d_y**2))
                    elif is_complex == 1:
                        tot_imag += h_y * (
                                -3*d_y**2*(a_0+a_2) +
                                 2*d_y   *(a_0-a_2)*h_y +
                                   h_y**2*(a_0+4*a_1+a_2)
                        ) / (3*(h_y**2-d_y**2))

        return [tot_real, tot_imag]

    elif method == 'separate':  # equidistant is assumed
        '''
        Instead of dealing with the summation over the coefficients, in this
        implementation the pressure function is interpolated by a cubic spline.
        Since the derivative of the Green's function is known analytically,
        we do not interpolate it in this implementation, however, this does
        need to be done for faster evaluation. Afterwards we calculate the integral
        over the subinterval by means of a regular sum (thus diving it into
        subsubintervals). Once again, we do not use a summation over coefficients,
        since this was easier to implement and produces approximately the same
        results.
        '''
        ### First we interpolate the pressure graph
        real_pressure = interp2d(meas_x, meas_y,
                            p(x, y)[0], kind='cubic')(plot_x, plot_y)
        imag_pressure = interp2d(meas_x, meas_y,
                            p(x, y)[1], kind='cubic')(plot_x, plot_y)
        ### Then for the greens function we use a finer grid
        real_greens = interp2d(plot_x, plot_y,
                        dG(plot_X, plot_Y, x_A, y_A)[0], kind='cubic')(plot_x, plot_y)
        imag_greens = interp2d(plot_x, plot_y,
                        dG(plot_X, plot_Y, x_A, y_A)[1], kind='cubic')(plot_x, plot_y)

        real = real_pressure*real_greens - imag_pressure*imag_greens
        imag = real_pressure*imag_greens + imag_pressure*real_greens
        res = [real, imag]

        return integral(plot_X, plot_Y, res, method='trapezeoid')

    elif method == 'altered_separate':  # non-equidistant is assumed
        '''
        See "separate" method
        '''
        ### First we interpolate the pressure graph
        points = np.array([x.flatten(), y.flatten()]).transpose()

        # We get an error if values to interpolate to not lie inside the convex
        # hull of points, thus the outer points are assumed to be equidistant
        # else, a nearest value approach is recommended.
        xmeas = x.copy()
        xmeas[1:-1,1:-1] = meas_X[1:-1,1:-1].copy()
        ymeas = y.copy()
        ymeas[1:-1,1:-1] = meas_Y[1:-1,1:-1].copy()

        interpol = np.array([xmeas.flatten(), ymeas.flatten()]).transpose()
        a = griddata(points, p(x, y)[0].flatten(), interpol, method='cubic')
        b = griddata(points, p(x, y)[1].flatten(), interpol, method='cubic')
        a = a.reshape((samples_y, samples_x))
        b = b.reshape((samples_y, samples_x))

        real_pressure = interp2d(meas_x, meas_y, a, kind='cubic')(plot_x, plot_y)
        imag_pressure = interp2d(meas_x, meas_y, b, kind='cubic')(plot_x, plot_y)
        ### Then for the greens function we use a finer grid
        real_greens = interp2d(plot_x, plot_y,
                dG(plot_X, plot_Y, x_A, y_A)[0], kind='cubic')(plot_x, plot_y)
        imag_greens = interp2d(plot_x, plot_y,
                dG(plot_X, plot_Y, x_A, y_A)[1], kind='cubic')(plot_x, plot_y)

        real = real_pressure*real_greens - imag_pressure*imag_greens
        imag = real_pressure*imag_greens + imag_pressure*real_greens
        res = [real, imag]

        return integral(plot_X, plot_Y, res, method='trapezeoid')


### defining variables
inter_x, inter_y = 50, 70
samples_x, samples_y = 29, 29
plot_samples_x, plot_samples_y = 4001, 4001
sources = [(10, 0, 2), (0, 15, 1.5), (1, -5, 1.7), (-13, 13, 2.3)]
k = 0.5  # omega/c, constant here
x_A, y_A, z_A = 0.9, -1.4, 10


### making measurement ticks
meas_x = np.linspace(-inter_x/2, inter_x/2, samples_x)
meas_y = np.linspace(-inter_y/2, inter_y/2, samples_y)
meas_X, meas_Y = np.meshgrid(meas_x, meas_y)
# meas_res = prod(meas_X, meas_Y, x_A, y_A)

## with noise
noise = 0.0005  # relative noise
np.random.seed(1)
meas_X2, meas_Y2 = meas_X.copy(), meas_Y.copy()
x_noise = np.random.normal(size=(samples_y, samples_x))*(inter_x*noise/samples_x)
y_noise = np.random.normal(size=(samples_y, samples_x))*(inter_y*noise/samples_y)
meas_X2 += x_noise
meas_Y2 += y_noise
meas_res = prod(meas_X2, meas_Y2, x_A, y_A)

### printing values
h_x = meas_x[1]-meas_x[0]
h_y = meas_y[1]-meas_y[0]
amount_noise_x = np.abs(x_noise).sum()/(samples_x*samples_y)
amount_noise_y = np.abs(y_noise).sum()/(samples_x*samples_y)
print(f'relative x noise of {amount_noise_x/h_x*100:.2g}%')
print(f'relative y noise of {amount_noise_y/h_y*100:.2g}%')

### making plot ticks
plot_x = np.linspace(-inter_x/2, inter_x/2, plot_samples_x)
plot_y = np.linspace(-inter_y/2, inter_y/2, plot_samples_y)
plot_X, plot_Y = np.meshgrid(plot_x, plot_y)
plot_res = prod(plot_X, plot_Y, x_A, y_A)

# calculate the answer for the integral on a finer grid
ans = integral(plot_X, plot_Y, plot_res)

# approximate the integral using the basic method
compare = 'trapezeoid'
val = integral(meas_X2, meas_Y2, meas_res, method=compare)
compare_rel_mse = np.sqrt(
    1/2*(abs(val[0]-ans[0])**2 + abs(val[1]-ans[1])**2) / (ans[0]**2 + ans[1]**2)
)
print(f'{compare.ljust(20)} method, with a MSE of: {compare_rel_mse:#.4e}')

# approximate the integral using various methods
for method in ['altered_trapezeoid', 'simpson',
               'altered_simpson', 'separate', 'altered_separate']:
    val = integral(meas_X2, meas_Y2, meas_res, method=method)
    rel_mse = np.sqrt(
        1/2*(abs(val[0]-ans[0])**2 + abs(val[1]-ans[1])**2) / (ans[0]**2+ans[1]**2)
    )
    print(f'{method.ljust(20)} method, with a MSE of: {rel_mse:#.4e} '
          f'this is a factor {compare_rel_mse/rel_mse:#.3g} better')
