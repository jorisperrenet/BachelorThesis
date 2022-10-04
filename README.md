# Methods for reducing error in approximations of the Rayleigh integral
This is the uncompiled version of my Bachelor Thesis at the TUDelft.
(The compiled version can be viewed here as the file `Approximating_Rayleigh_Integrals_Joris_Perrenet.pdf` and in the [repository](http://resolver.tudelft.nl/uuid:c48fa27f-d91b-4c07-a74b-6ac6762bc095))

Also, a huge thanks to [VimTeX](https://github.com/lervag/vimtex) for making it so easy to compile my tex files and interacting with [zathura](https://github.com/pwmt/zathura) enabling hot reload. It allowed me to edit my thesis in [NeoVim](https://github.com/neovim/neovim) whilst only needing to save the file in order to see the changes in my pdf-viewer.

Some notes:
1.  To compile my thesis I used the compiler `latexmk` (it completely compiles without any errors and warnings) after installing [`texlive-most`](https://archlinux.org/groups/x86_64/texlive-most/) on Arch Linux.
2.  The `TUD-report2020.cls` and `codestyle.sty` in the `report` folder can also be used as templates for other people at the TUDelft (these are not the standard files, I edited them a lot).
3.  Not all used packages are in `Approximating_Rayleigh_Integrals_Joris_Perrenet.pdf`, some are in `TUD-report2020.cls` as it was needed to ensure that everything was loaded in the correct order.
4.  The lemmas and theorems used in appendix B of my thesis are also listed explicitely in the `lemmas` folder.
5.  All the used `.tikz` files are in `report/pictures`.
6.  To export `.svg` files I used the following command:
```{bash}
inkscape -D -z --file=<FILENAME.svg> --export-pdf=<FILENAME.pdf> --export-latex
```
I highly recommend using `.svg` (and `.tikz`) files as they are vector files and extremely lightweight, it allowed my 73-paged thesis to be ~1.13 MB.
