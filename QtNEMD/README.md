# QtNEMD - A friendly python front-end for non-equilibrium molecular dynamics
QtNEMD is a python application which provides a friendly graphical front-end in Qt for molecular 
dynamics calculations using software for developed by the Bernhardt group at the University of
Queensland.

This project is currently a work in progress, bug reports and feature requests welcome.

# Installation
## Prerequisites
You'll need the following software packages installed on your computer:
  1) Python 3
  2) The `pip` python package manager (usually bundled with your python distribution)
  3) `make` (the software has been tested on Gnu Make, but other dialects may work)
  4) The `gfortran` compiler (although the `Makefile` can be edited to use other compilers)

## Installation with `pip`
This software depends on the following python modules:
  1) `numpy`
  2) `f2py` Fortran compatibility layer (part of numpy)
  3) `PyQt5` graphical user interface library
  4) `pyqtgraph`

These packages can be automatically installed by `pip` via the following command:

```
pip install -e .
```

This will run the `setup.py` script, which installs the correct packages in a new folder in the root 
directory of the repository. `setup.py` also calls `make` to build the Fortran backend and GUI assets.

## Manual installation with `make`
If you already have the required python packages installed, you can simply run `make` to build the
Fortran backend and GUI assets. It is also possible to build these components separately by running
`make driver` for Fortran and `make gui` for GUI assets. Note that graphical assets are automatically
build by the `pyuic5` tool, which should be included in `PyQt` if you installed via `pip`, but you may 
need to install manually it if you used your OS's package manager. Finally, `make clean` removes the
compiled Fortran files and GUI assets if you need to do a clean build.
