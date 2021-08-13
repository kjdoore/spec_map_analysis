## General Description

**Spec_map_analysis** is a spectral map analysis package. The package contains modules that allow for the fitting of emission lines in spectral maps, resampling and convolution of images, and derivation of emission line ratios and subsequently electron densities. 

 - The emission line fitting module can fit a single emission line with a Gaussian profile and optionally a polynomial continuum component. The spectral fitting is performed independently for each spatial pixel using a non-linear least squares method.

- The resampling and convolution module resamples one image to another image's pixel resolution, and then convolves either image to the specified spatial resolution. The convolution allows for the propagation of uncertainties, assuming they are correlated, using the method presented in [Klein (2021)](https://iopscience.iop.org/article/10.3847/2515-5172/abe8df). 

- The line ratio and electron density module allow for the derivation of specified line ratios and subsequently, electron density maps using atomic data from [PyNeb](http://research.iac.es/proyecto/PyNeb//). Aperture photometry can also be performed on the line intensity maps, and then used to calculate line ratios and electron densities for aperture regions.


## Requirements

**Spec_map_analysis** has the following requirements:

 - [python](https://www.python.org) 3.6 or later
 - [numpy](https://numpy.org) 1.20 or later
 - [astropy](https://www.astropy.org) 4.2 or later
 - [scipy](https://www.scipy.org) 1.6 or later
 - [reproject](https://reproject.readthedocs.io/en/stable/#) 0.7 or later
 - [PyNeb](http://research.iac.es/proyecto/PyNeb//) 1.1 or later

**Spec_map_analysis** also has optionally dependencies on:

 - [photutils](https://photutils.readthedocs.io/en/stable/index.html) 1.1 or later: To derive aperture photometry.
 - [matplotlib](https://matplotlib.org) 3.3 or later: To plot images.
 - [emcee](https://emcee.readthedocs.io/en/stable/) 3.0 or later: To fit spectral maps using an MCMC algorithm.


## Installation

**Spec_map_analysis** can be downloaded, installed and upgraded using [pip](https://pip.pypa.io/en/latest/).

To install only the required packages, run:
```
pip install spec_map_analysis
```
To install all package dependencies, run:
```
pip install spec_map_analysis[all_dep]
```
To install directly from the [GitHub repository](https://github.com/kjdoore/spec_map_analysis), run:
```
pip install git+https://github.com/kjdoore/spec_map_analysis.git
```


## Manual

The manual can be found [here](https://github.com/kjdoore/spec_map_analysis/tree/master/docs).


## Reporting Issues

If you happen to find a bug in spec_map_analysis, please report it by creating a new issue on the [**spec_map_analysis** GitHub issue tracker](https://github.com/kjdoore/spec_map_analysis/issues). Please include an example of the issue that allows for its reproduction, so the problem can be fixed.