# Resampling and Convolution ([spec_map_analysis.reproject_convolve](https://github.com/kjdoore/spec_map_analysis/blob/master/docs/user_docs/reproject_convolve.md))

## Introduction

In **spec_map_analysis**, the **reproject_convolve** sub-package is used to resample (reproject) and convolve two images and their corresponding uncertainty images to the same pixel scale and spatial resolution. The convolution of the uncertainties utilizes the method presented in [Klein (2021)](https://iopscience.iop.org/article/10.3847/2515-5172/abe8df) that allows for proper propagation of uncertainties, assuming correlation between pixels. The result of the convolution and resampling is two FITS files each containing the convolved and/or resampled image extensions of the data and corresponding uncertainty.


## Map Resampling and Convolution 

To resample and convolve two images, a single call to the [**reproject_convolve()**](https://github.com/kjdoore/spec_map_analysis/blob/master/spec_map_analysis/reproject_convolve/reproject_convolve.py) function is needed. The required inputs include the names of the two output files, the names of the two input files, and the spatial resolution of each image. The first input image must be the higher spatial resolution image that will be convolved to the lower spatial resolution. The second input image gets resampled by default, but that can be changed with the keyword `file_to_reproject`. A value of 1 indicates the first input image should be resampled, and a value of 2 (default) indicates the second. 

To give an example, we can use an [OIII] 52&mu;m line intensity map from SOFIA/FIFI-LS and an [OIII] 88&mu;m line intensity map from *Herschel*/PACS of the same source derived with the [**spec_map_analysis.spectra_fitting.line_map()**](https://github.com/kjdoore/spec_map_analysis/blob/master/docs/user_docs/spectra_fitting.md) function. The spatial resolution of these maps would be 6.2'' and 9'' in terms of FWHM, respectively. To convolve the FIFI-LS image to the PACS spatial resolution and resample the PACS image to the FIFI-LS pixel, we can do:

```python
from spec_map_analysis.reproject_convolve import reproject_convolve
import numpy as np

input_file1 = 'line_intensity_03_52.fits'
input_file2 = 'line_intensity_03_88.fits'

output_file1 = 'line_intensity_03_52_convolved.fits'
output_file2 = 'line_intensity_03_88_resampled.fits'

# Divide the FWHM by 2*sqrt(2*ln(2)) to convert FWHM to sigma
reproject_convolve(output_file1, output_file2, input_file1, input_file2, 6.2 / (2 * np.sqrt(2 * np.log(2))),
                   9 / (2 * np.sqrt(2 * np.log(2))))
```

If we instead wanted to both convolve and resample the FIFI-LS image, we could simply set `file_to_reproject=1`. Also, to make things simpler for the user, the keyword `fwhm` can be set to `True`, which allows the spatial resolution to be input as a FWHM rather than a standard deviation.

**NOTE**: The data and uncertainty extensions for the input files are not required inputs, rather they are optional keywords that can be set. By default, the extensions are set to be the line intensity and uncertainty extensions produced by the [**spec_map_analysis.spectra_fitting.line_map()**](https://github.com/kjdoore/spec_map_analysis/blob/master/docs/user_docs/spectra_fitting.md) function. If other data and uncertainty extensions are used, the `data_extension` and `unc_extension` keywords can be set to the data and uncertainty extensions, respectively. Note that the extensions for the input files must be the same.

### Pixel Unit Handling

When resampling an image, it is important to properly handle the units. For example, a flux in **Jy** in a given pixel will either be per pixel (**Jy/pixel**) or per an angular unit (e.g., **Jy/deg**). Depending on whether the unit is per pixels or per angular unit, the values in each pixel may need to be rescaled. Per pixel units will need to be rescaled to account for the change in pixel size when resampling, while per angular unit will not. Since input data can have either unit type, [**reproject_convolve()**](https://github.com/kjdoore/spec_map_analysis/blob/master/spec_map_analysis/reproject_convolve/reproject_convolve.py) has the keyword `pixel_unit` to handle any rescaling of the data if needed. The keyword values are `pixel` (default) or `angular`, which mean the input units are per pixel or per angular unit, respectively.


### Resampling Method

Different algorithms exist for resampling an image. To use a different algorithm other than the default linear interpolation algorithm, the keyword `reprojection_method` can be set. The options are:

- `'adaptive'`
- `'exact'`
- `'interp'` (default)

The value of the `reprojection_method` keyword determines what algorithm is used. Setting the keyword utilizes a different function in the [reproject](https://reproject.readthedocs.io/en/stable/#) package with the corresponding name. Details on each function can be found at <https://reproject.readthedocs.io/en/stable/celestial.html>