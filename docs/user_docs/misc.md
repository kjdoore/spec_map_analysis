# Miscellaneous ([spec_map_analysis.misc](https://github.com/kjdoore/spec_map_analysis/blob/master/docs/user_docs/misc.md))

## Introduction

In **spec_map_analysis**, the **misc** sub-package contains miscellaneous functions that may be useful when deriving line intensity, line ratio, or electron density maps. 


## Effective Exposure Times

When fitting spectral maps with exposure times that vary across the image, it may be necessary to calculate the effective exposure time of the fit emission line. This can easily be done with the [**effective_exp_time()**](https://github.com/kjdoore/spec_map_analysis/blob/master/spec_map_analysis/misc/effective_exp_time.py) function. The function takes the original spectral cube file and the resulting fit from the [**spec_map_analysis.spectra_fitting.line_map()**](https://github.com/kjdoore/spec_map_analysis/blob/master/docs/user_docs/spectra_fitting.md) function to compute the effective exposure time within the FWHM of the fit line.

For example, we can determine the effective exposure time for the fit a spectral cube using:

```python
from spec_map_analysis.misc import effective_exp_time

input_file = 'spectral_cube.fits'
fit_file = 'line_intensity.fits'
out_file = 'effective_exp_time.fits'
# Extension names of the input file
extensions = ['FLUX', 'ERROR', 'WAVELENGTH', 'EXPOSURE_MAP']

effective_exp_time(out_file, input_file, extensions, fit_file)
```

The resulting output file will contain an image extension with the effective exposure time of the fit lines.


## Cubik Uncertainties

The **spec_map_analysis** package was originally designed to fit spectra and derive electron density maps for SOFIA/FIFI-LS spectral maps. However, spectral maps output by the FIFI-LS pipeline tend to have underestimated uncertainties. To get better estimates for the FIFI-LS uncertainties, the code [cubik](https://github.com/darioflute/fifipy), by [Dario Fadda](https://github.com/darioflute/), can be used to generate better estimates using a spectrally and spatially weighted average of the raw data points.

The [**cubik_uncertainty()**](https://github.com/kjdoore/spec_map_analysis/blob/master/spec_map_analysis/misc/cubik_uncertainty.py) function was made to simply replace the uncertainties in the pipeline file with those produced by [cubik](https://github.com/darioflute/fifipy). This can be done with the simple example:

```python
from spec_map_analysis.misc import cubik_uncertainty

pipeline_file = 'spectral_cube.fits'
cubik_file = 'WXY_cubik.fits'
output_file = 'spectral_cube_cubik.fits'

# Extensions for the pipeline and cubik errors
err_extensions = ['ERROR', 'ERROR']

cubik_uncertainty(output_file, pipeline_file, cubik_file, err_extensions)
```






