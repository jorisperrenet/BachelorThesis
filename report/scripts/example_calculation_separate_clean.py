'''
In this script we approximate the integral of f(x)*g(x) over the
interval [-x_L, x_R]. The used method is described in chapter 3.1
and the used alteration to the method is given in chapter 3.3.
'''
import numpy as np
import numpy.polynomial.polynomial as np_pol
from scipy.interpolate import interp1d


### We first define both functions
def f(x):
    return np.cos(np.sqrt((x-1)**2/4+1))
    # return 1/(1+(abs(x)**2)/10)

def g(x):
    return np.sin(np.sqrt((x-2)**2/4+1))
    # return 1/(1+(abs(x+1)**2)/10)


### Here we plot the difference between separate and combined spline interpolation
x_L, x_R = -6, 6
n_f, n_g = 3, 5
num_samples = 9

### Defining the samples points and the distances/inner_widths
### of the subintervals, i.e. x_1-x_0
ell = num_samples-1
samples = np.linspace(x_L, x_R, num_samples)
inner_width = (x_R-x_L)/(ell)

### Interpolating the not known polynomial with (cubic) splines at the
### sample points
splines = interp1d(samples, f(samples), kind='cubic')

### Calculating the coefficients of interpolation, i.e. matrices B and C
# Note: The function interp1d from scipy does not give the coefficients of the
#       cubic splines therefore we use numpy for a polynomial fit on the cubic
#       splines on all the subintervals from the fit we can extrapolate the
#       coefficients of the polynomials. However, this method is not recommended
#       as it requires a subdivision of the interval, instead implement a spline
#       interpolation that does return the interpolating coefficients.
b, c = [], []
for idx, i in enumerate(samples[:-1]):
    # with shifted intervals
    interval = np.linspace(i, i+inner_width, 1000)
    mapped_interval = np.linspace(0, inner_width, 1000)
    pol_f = np_pol.Polynomial.fit(mapped_interval, splines(interval), n_f)
    pol_g = np_pol.Polynomial.fit(mapped_interval, g(interval), n_g)

    b.append(list(pol_f.convert().coef))
    c.append(list(pol_g.convert().coef))

### Calculating and storing the powers of x, note that the intervals were shifted
powers_x = [
        [(inner_width)**(jk+1) / (jk+1) for jk in range(n_f+n_g+2)]
        for i in range(len(c))
]

### Calculating and storing all values of the matrix D
# Note: We do not need to store the values of the matrix D if we only have
#       one sensor because we then do not need to shift the matrix B to calculate
#       its frobenius product with matrix D, however for clarity this is
#       implemented nonetheless.
d = [[0 for _ in range(n_f+1)] for _ in range(ell)]
for i in range(ell):
    for j in range(n_f+1):
        d[i][j] = sum(c[i][k] * powers_x[i][j+k] for k in range(n_g+1))

### Calculating the final integral
som = 0
for i in range(ell):
    for j in range(n_f+1):
        som += b[i][j] * d[i][j]

print(f'Our answer for the integral is {som}')
