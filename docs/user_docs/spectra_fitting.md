# Spectral Map Fitting ([spec_map_analysis.spectra_fitting](https://github.com/kjdoore/spec_map_analysis/blob/master/docs/user_docs/spectra_fitting.md))

## Introduction

In **spec_map_analysis**, the **spectra_fitting** sub-package is the emission line fitting procedures used to fit a spectral cube containing a single emission line. The emission line is fit with a Gaussian profile and can optionally include a polynomial continuum component up to degree two that is fit simultaneously with the Gaussian profile. The spectral fitting is performed independently for each spatial pixel using a non-linear least squares method. Ranges for the fitting parameters can be limited to user specified values, automatically limited based off of an initial fit's high signal-to-noise parameter range, or both. The result is a FITS (Flexible Image Transport System) file containing the image extensions of the fit line intensity and corresponding uncertainty, continuum emission, velocity offset, full width at half maximum (FWHM) of the line, and fit parameters and their covariance matrix. 


## Fitting Spectral Maps

Fitting a spectral map with default settings only requires a simple call of the [**line_map()**](https://github.com/kjdoore/spec_map_analysis/blob/master/spec_map_analysis/spectra_fitting/line_map.py) function. The only required inputs are the name of output FITS file; the input FITS file; a list of extensions for the input file that contain the spectral cube, its uncertainty, and wavelengths of the spectra; and the peak rest-frame wavelength of the emission line.

**NOTE**: The input data must be in units of **Jy** for the fluxes and **&mu;m** for the wavelengths. No options for unit conversion are provided.

For example, we can fit a spectral cube containing the IR emission line [OIII] 52&mu;m using:
```python
from spec_map_analysis.spectra_fitting import line_map

input_file = 'spectral_cube.fits'
out_file = 'line_intensity.fits'

# Extension names of the input file
extensions = ['FLUX', 'ERROR', 'WAVELENGTH']

# Peak rest-frame wavelength of the [OIII] emission line
center_wavelength = 51.8145

line_map(out_file, input_file, extensions, center_wavelength)
```

However, this assumes that the wavelengths in the wavelength extension are in the galaxy's rest-frame. If this is not the case, the velocity of the galaxy for a chosen standard of rest can be provided with the `velocity` keyword, which must be in units of **km/s**.

The result from the fitting is saved to the given `output_file` as a FITS file that contains the following image extensions: the fit line intensity, the fit line intensity uncertainty, the continuum emission at the peak of the emission line, the velocity offset of the peak of the fit line from the expected central peak, the FWHM of the line, the best fit parameters, and the best fit parameters covariance matrix. 

### Fitting with Continuum

Fitting the emission line simultaneously with continuum emission simply requires the addition of the `nterms` keyword in the call of [**line_map()**](https://github.com/kjdoore/spec_map_analysis/blob/master/spec_map_analysis/spectra_fitting/line_map.py). The values of `nterms` must be integer values from 3 to 6, where each value has the following continuum component:

- `nterms=3` - Fits the emission line **without** a continuum component (Default)

- `nterms=4` - Fits the emission line with a **constant or zero degree polynomial** continuum component

- `nterms=5` - Fits the emission line with a **linear or first degree polynomial** continuum component
  
- `nterms=6` - Fits the emission line with a **quadratic or second degree polynomial** continuum component

The name `nterms` comes from the number of terms in the function used to fit the data. The function has the form

f(x) = a<sub>0</sub> exp(-0.5 ((x-a<sub>1</sub>)/a<sub>2</sub>)<sup>2</sup>) + a<sub>3</sub> + a<sub>4</sub>x +a<sub>5</sub>x<sup>2</sup>,

where a<sub>0</sub> is the peak height of the fit Gaussian, a<sub>1</sub> is the peak offset, a<sub>2</sub> is the sigma of the Gaussian, and a<sub>3</sub>, a<sub>4</sub>, and a<sub>5</sub> are the continuum polynomial terms. Therefore, when `nterms=3`, only the parameters of the Gaussian, a<sub>0</sub>, a<sub>1</sub>, and a<sub>2</sub>, are fit, and when `nterms=4`, the constant continuum term a<sub>3</sub> is also fit, and so on. 

### Fitting with Specified or Revised Parameter Ranges

The [**line_map()**](https://github.com/kjdoore/spec_map_analysis/blob/master/spec_map_analysis/spectra_fitting/line_map.py) function has the capability for the allowed ranges of fitting parameters to be user specified to prevent unrealistic values. Alternatively, parameter ranges can be automatically restricted based off of the parameter range of high signal-to-noise pixels from an initial fit to the data. Both the user specified and automatic methods can be used simultaneously, which will only restrict the user specified parameter range further.

To manually set the parameter ranges, the keywords `lower_bounds` and `upper_bounds` must be set in tandem to `numpy` arrays with lengths of `nterms`. Each value in the arrays correspond to the parameters given by the function defined by `nterms`. Therefore, the first value must be the min/max height of the Gaussian in **Jy**; the second is the min/max peak offset in **km/s**; the third is the min/max FWHM of the Gaussian in **km/s**, which is converted to sigma; and the min/max polynomial terms in **Jy** depending on the value of `nterms`. For a specific parameter to have an unlimited range, the value in `lower_bounds` or `upper_bounds` must be set to negative or positive infinity, respectively, using the `numpy` infinity values (`numpy.inf`).

For example to fit for an emission line limiting the peak offset and FWHM of the Gaussian to be within +/-100 km/s and under 350 km/s, respectively, without limiting the emission peak height and without continuum, we can run:

```python
from spec_map_analysis.spectra_fitting import line_map
import numpy as np

input_file = 'spectral_cube.fits'
out_file = 'line_intensity.fits'

# Extension names of the input file
extensions = ['FLUX', 'ERROR', 'WAVELENGTH']
# Peak rest-frame wavelength of the [OIII] emission line
center_wavelength = 51.8145

# Parameter range limits
upper_bounds = np.array([np.inf, 100, 350])
lower_bounds = np.array([0, -100, 0])

line_map(out_file, input_file, extensions, center_wavelength, 
         lower_bounds=lower_bounds, upper_bounds=upper_bounds)
```

To automatically restricted the parameter ranges, `revise_bounds` must be set to `True`. If so, an initial fit is performed on the data. Then, the range of fit parameters from pixels that have SNR values above the value specified in `snr_limit` are used as the new parameter ranges, and the data is fit again for the final result.

For example, to automatically restrict the parameter ranges to that of the pixels with SNR > 3,  we run:

```python
from spec_map_analysis.spectra_fitting import line_map

input_file = 'spectral_cube.fits'
out_file = 'line_intensity.fits'

# Extension names of the input file
extensions = ['FLUX', 'ERROR', 'WAVELENGTH']

# Peak rest-frame wavelength of the [OIII] emission line
center_wavelength = 51.8145

line_map(out_file, input_file, extensions, center_wavelength, 
         revise_bounds=True, snr_limit=3)
```


### Fitting with the MCMC Algorithm

The **spectra_fitting** sub-package also allows for the fitting of the emission lines using the Affine Invariant Markov chain Monte Carlo (MCMC) Ensemble sampler as implemented in the [emcee package](https://emcee.readthedocs.io/en/stable/#). To fit the data with the MCMC algorithm, simply call the [**line_map()**](https://github.com/kjdoore/spec_map_analysis/blob/master/spec_map_analysis/spectra_fitting/line_map.py) function as shown in any of the examples above, but now with the `mcmc` keyword set as `True`. The keywords specific to the MCMC algorithm (`nwalkers`, `nsteps`, and `discard`) can optionally be set to change the number of walkers in the MCMC ensemble, the number of steps taken in the MCMC chain, and the number of burn-in phase steps that should be discarded, respectively.

For example, to fit a spectral map of [OIII] 52 &mu;m without continuum emission and unrevised parameter ranges, using the MCMC algorithm, run:

```python
from spec_map_analysis.spectra_fitting import line_map

input_file = 'spectral_cube.fits'
out_file = 'line_intensity.fits'

# Extension names of the input file
extensions = ['FLUX', 'ERROR', 'WAVELENGTH']

# Peak rest-frame wavelength of the [OIII] emission line
center_wavelength = 51.8145

line_map(out_file, input_file, extensions, center_wavelength, 
         mcmc=True)
```

The result from the MCMC fitting is the same as the standard output, except the best fit parameters and covariance matrix extensions are replaced with the MCMC chain for each parameter.

### Note

See the description of [**line_map()**](https://github.com/kjdoore/spec_map_analysis/blob/master/spec_map_analysis/spectra_fitting/line_map.py) and [**line_fitting()**](https://github.com/kjdoore/spec_map_analysis/blob/master/spec_map_analysis/spectra_fitting/line_fitting.py) for a full list and detailed description of all optional parameters not discussed here.