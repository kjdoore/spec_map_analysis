# Line Intensity Upper Limit Estimates ([spec_map_analysis.upper_limits](https://github.com/kjdoore/spec_map_analysis/blob/master/docs/user_docs/upper_limits.md))

## Introduction

In **spec_map_analysis**, the **upper_limits** sub-package is a method for deriving a first order estimate of the upper limit line intensity of a given spectral map. The method consists of randomly inserting simulated sources into a spectral map and determining at what peak height and subsequent line intensity are the simulated lines detected above 1&sigma;. The result consists of an array of these detection limit line intensities that can be used to determine a median and error range on the upper limit.


## Estimating Upper Limits

To generate the array of upper limit estimates, only a single call to the [**upper_limit_estimate()**](https://github.com/kjdoore/spec_map_analysis/blob/master/spec_map_analysis/upper_limits/upper_limit_estimate.py) function is needed. The function iteratively inserts a source at a random location. It then increases the height of the simulated source at that location until the line is detected by the fitting algorithm, or the maximum allowed height of the simulated source is reached. The width and velocity offset of the simulated source are user supplied, and the number of inserted sources can optionally be set using the `nsim` keyword.

**NOTE**: The units for the maximum peak height, velocities, and FWHM must be in **Jy**, **km/s**, and **km/s**, respectively. No options for unit conversion are provided.

For example, we can estimate the upper limit array of a spectral cube containing the IR emission line [OIII] 52&mu;m with 50 simulated sources using:
```python
from spec_map_analysis.upper_limits import upper_limit_estimate as ule
input_file = 'spectral_cube.fits'
# Extension names of the input file
extensions = ['FLUX', 'ERROR', 'WAVELENGTH']
# Peak rest-frame wavelength of the [OIII] emission line
center_wavelength = 51.8145
# Velocity of the galaxy to the chosen standard of rest in km/s
velocity = 0
# Maximum simulated source peak height in Jy
max_height = 1
# Simulated source velocity offset and FWHM in km/s
peak_velocity = 100
fwhm = 350

upper_limits = ule(input_file, extensions, center_wavelength, velocity, 
                   max_height, peak_velocity, fwhm, nsim=50)
```


### Excluding a Masked Region

The [**upper_limit_estimate()**](https://github.com/kjdoore/spec_map_analysis/blob/master/spec_map_analysis/upper_limits/upper_limit_estimate.py) function can also mask a region within the map when inserting simulated sources. This mask prevents simulated sources from being placed in the mask region. Masking locations of known sources will allow for a better estimate on the upper limit.

This can be done for the example above using:
```python
from spec_map_analysis.upper_limits import upper_limit_estimate as ule
from astropy.io import fits
input_file = 'spectral_cube.fits'
# Extension names of the input file
extensions = ['FLUX', 'ERROR', 'WAVELENGTH']
# Peak rest-frame wavelength of the [OIII] emission line
center_wavelength = 51.8145
# Velocity of the galaxy to the chosen standard of rest in km/s
velocity = 0
# Maximum simulated source peak height in Jy
max_height = 1
# Simulated source velocity offset and FWHM in km/s
peak_velocity = 100
fwhm = 350
# Read in data from mask file
mask = fits.open('mask_region.fits').data

upper_limits = ule(input_file, extensions, center_wavelength, velocity, 
                   max_height, peak_velocity, fwhm, mask=mask)
```


### Note

See the description of [**upper_limit_estimate()**](https://github.com/kjdoore/spec_map_analysis/blob/master/spec_map_analysis/upper_limits/upper_limit_estimate.py) for a full list and detailed description of all optional parameters not discussed here, which include the `nterms`, `lower_bounds`, and `upper_bounds` parameters discussed in [**spec_map_analysis.spectra_fitting**](https://github.com/kjdoore/spec_map_analysis/blob/master/docs/user_docs/spectra_fitting.md).