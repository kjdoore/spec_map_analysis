# Plotting Images ([spec_map_analysis.plots](https://github.com/kjdoore/spec_map_analysis/blob/master/docs/user_docs/plots.md))

## Introduction

In **spec_map_analysis**, the **plots** sub-package is used to plot spectra and images. Basic spectra from the spectral map at a given location can be plotted, and additionally the fit to the spectrum can be included with the resulting residual. The function used to plot images is highly versatile and allows for any number of images to be plotted in the same graphic. The only limitation is that the images must be at a similar sky location as all images will share the same axes range. The sub-package utilizes [matplotlib](https://matplotlib.org) and saves the resulting graphics to a specified file. 


## Spectra Plots

The [**spectrum_plot()**](https://github.com/kjdoore/spec_map_analysis/blob/master/spec_map_analysis/plots/spectrum_plot.py) function is specifically designed to plot multiple spectra from a spectral map at different specified pixel locations. Locations of the spectra can be specified in either pixel coordinates or sky coordinates (e.g., Right Ascension and Declination). 

Since the locations of the spectra must be specified, the file name to save the graphics will be an initial name that will be appendend to specify the location of the spectrum. For example, we can plot the spectra at pixel locations (2, 5) and (10, 3) using:

```python
from spec_map_analysis.plots import spectrum_plot

input_file = 'spectral_cube.fits'
out_file = 'spectra'
# Extension names of the input file need for plotting
extensions = ['FLUX', 'ERROR', 'WAVELENGTH']
# Pixel locations of the spectra to plot
locations = [(2, 5), (10, 3)]

spectrum_plot(out_file, input_file, extensions, locations)
```

The resulting saved graphics will have the following names: `spectra_pixel_x2_y5.png` and `spectra_pixel_x10_y3.png`.

If sky coordinates are given instead of pixel coordinates, then the `location_unit` keyword must be set to `'sky'` to signal that the locations are given in sky coordinates. Coordinates given this way must have the RA and Dec in seperate strings and in a format understood by [astropy.coordinates.SkyCoord](https://docs.astropy.org/en/stable/api/astropy.coordinates.SkyCoord.html). Additional keywords (kwargs) are allowed in [**spectrum_plot()**](https://github.com/kjdoore/spec_map_analysis/blob/master/spec_map_analysis/plots/spectrum_plot.py) and are supplied to [astropy.coordinates.SkyCoord](https://docs.astropy.org/en/stable/api/astropy.coordinates.SkyCoord.html), if for example, the specified sky coordinates are in a frame other than ICRS. Note that the locations supplied to [**spectrum_plot()**](https://github.com/kjdoore/spec_map_analysis/blob/master/spec_map_analysis/plots/spectrum_plot.py) must be a list to two element tuples, regardless of the number of locations and coordinate type.

A secondary y-axis can be added to the plotted spectrum showing the velocity offset from the rest-frame in **km/s**. To add the secondary axis, the `center_wavelength` keyword must be set to the expected center wavelength of the spectral line as if observed in the rest-frame. The units of `center_wavelength` must be in **&mu;m** or match the wavelength unit specified in the `units` keyword. Further, the `velocity` keyword should be specified as the velocity of the object in **km/s** to the chosen standard of rest if the velocity axis is to be in reference to the object's rest-frame.

To include the fit to the spectrum and the resulting residual, the `fit_file` keyword will need to contain the name of the fit file produced by the [**spec_map_analysis.spectra_fitting.line_map()**](https://github.com/kjdoore/spec_map_analysis/blob/master/docs/user_docs/spectra_fitting.md) function. To do this, we can update the above example to:

```python
from spec_map_analysis.plots import spectrum_plot

input_file = 'spectral_cube.fits'
out_file = 'spectra'
fit_file = 'line_intensity.fits'

# Extension names of the input file need for plotting
extensions = ['FLUX', 'ERROR', 'WAVELENGTH']
# Pixel locations of the spectra to plot
locations = [(2, 5), (10, 3)]

spectrum_plot(out_file, input_file, extensions, locations, fit_file=fit_file)
```

The resulting graphics will now contain a plot of the spectrum with the fit to the data overlayed and a plot of the resulting residual below it sharing the same x-axis.

By default, all graphics will be saved as .png files. To change the file type, the `out_file_type` keyword can be set to a three letter string, specifying the file type. File types must be supported by the [savefig](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.savefig.html) function in [matplotlib](https://matplotlib.org).


## Map Plots

The [**map_plot()**](https://github.com/kjdoore/spec_map_analysis/blob/master/spec_map_analysis/plots/map_plot.py) function is designed to be highly versatile and allow for any number of images to be plotted in the same graphic. The function takes an input file and image extension, plots the image, and saves the graphic to the specified output file. The function relies heavily on the flexibility of `kwargs` to allow for an arbitrary number of keyword arguments.


### Single Images

To plot a single image in the most simple example, we can use:

```python
from spec_map_analysis.plots import map_plot

input_file = 'line_intensity.fits'
out_file = 'intensity_map.png'
# Extension name of the image
extension = 'LINE_INTENSITY'

map_plot(out_file, input_file, extension)
```

This results in the line intensity image extension being plotted with a colorbar and saved to the output file. The resulting figure size may not be scaled properly, resulting in a crowded graphic. To adjust the figure size, the keywords `figure_height_scale` and `figure_width_scale` can be set to scalar values to scale the height and width of the figure as needed (default values are 1). For example, setting both keywords to a value of 2 will double the figure size in both height and width.

Additionally, the colormap of the plotted image can be set with the `cmap` keyword. The keyword can be any predefined colormap within [matplotlib](https://matplotlib.org) (a full list can be found [here](https://matplotlib.org/stable/tutorials/colors/colormaps.html)). The default value is the standard [matplotlib](https://matplotlib.org) default of `'viridis'`.

#### Adding Labels

The [**map_plot()**](https://github.com/kjdoore/spec_map_analysis/blob/master/spec_map_analysis/plots/map_plot.py) function will automatically read the FITS file header of the required input image extension. The image axis will then be set as the WCS axis as defined by the header with labels of RA and Dec. These labels are not changeable, but additional labels can be added for the colorbar and title of the image.

In the single image case, the keyword `cbar_label0` can be set to a string that labels the colorbar, and `title0` can be set to a string for the title of the image. These keywords, and all those ending with `0` below, are not specified in the list of keywords in the [**map_plot()**](https://github.com/kjdoore/spec_map_analysis/blob/master/spec_map_analysis/plots/map_plot.py) function, but instead will be passed to `kwargs`. The reasoning for this will become clear further below, when we want to plot multiple images.

Expanding the previous example to include these labels and a change in the colormap, we can use:

```python
from spec_map_analysis.plots import map_plot

input_file = 'line_intensity.fits'
out_file = 'intensity_map.png'
# Extension name of the image
extension = 'LINE_INTENSITY'

map_plot(out_file, input_file, extension, cmap='hsv',
         title0='Line Intensity Map', cbar_label0='W s$^{-2}$')
```

#### Setting Axis and Colorbar Limits

The axis range of the image, as stated above, is automatically determined by the WCS axis as defined by the FITS file header of the required input image extension. To manually adjust the axis limits, the keywords `ralim` and `declim` can be set to two element lists of strings containing the upper and lower Right Ascension and Declination limits, respectively. The string values must be clearly discernible by the [astropy.coordinates.Angle](https://docs.astropy.org/en/stable/api/astropy.coordinates.Angle.html#astropy.coordinates.Angle) function.

If the colorbar and data range that the colormap covers need to be limited to a set range, the `range0` keyword can be set. This keyword is a two element list giving the lower and upper limit of the data that is to be plotted. Data below and above the range will be plotted as the extrema colors of the colorbar.

#### Aperture Patches

Circular aperture patches can easily be overlayed on the plotted image in [**map_plot()**](https://github.com/kjdoore/spec_map_analysis/blob/master/spec_map_analysis/plots/map_plot.py). To overlay one or more aperture patches, the keyword `aperture0` can be set to a list of two element list of strings each containing the Right Ascension and Declination location to place the circular apertures. Like the `ralim` and `declim` keywords, the RA and Dec locations must be clearly discernible by the [astropy.coordinates.Angle](https://docs.astropy.org/en/stable/api/astropy.coordinates.Angle.html#astropy.coordinates.Angle) function. To specify the radius of each aperture, the `aperture_radius0` keyword will need to set as a list of scalars containing the corresponding radius in **arcseconds** of each aperture patch. Finally, the color of all apertures will be the same and can be specified as any acceptable color string to [matplotlib](https://matplotlib.org) using the `aperture_color0` keyword. If not specified, the apertures' colors will be red.

For example, we can add three different apertures to an image of NGC 1569. Two apertures with radii of 6'' and 12'' at the central peak and another with a radius of 6'' at the secondary emission peak. Here, we will also limit the axis range using:

```python
from spec_map_analysis.plots import map_plot

input_file = 'NGC1569_line_intensity.fits'
out_file = 'intensity_map_w_apertures.png'
# Extension name of the image
extension = 'LINE_INTENSITY'

map_plot(out_file, input_file, extension, title0='Line Intensity Map', cbar_label0='W s$^{-2}$',
         ralim=['4h30m56.5s', '4h30m42.5s'], declim=['64:50:17 degrees', '64:51:38 degrees'], 
         aperture0=[['4h30m47.0456s', '64d51m00.506s'], ['4h30m47.0456s', '64d51m00.506s'],
                    ['4h30m52.1933s', '64d50m46.827s']], aperture_radius0=[12 ,6, 6])
```

#### Contours

Contours of any image can also be easily overlayed on the plotted image in [**map_plot()**](https://github.com/kjdoore/spec_map_analysis/blob/master/spec_map_analysis/plots/map_plot.py).  To overlay contours, the keyword `contour0` can be set to a list of strings with each string containing a name of a FITS file whose contours are to be used. Also, the keyword `contour_extension0` will need to be set to a list containing the name or index of the extensions of each `contour0` file string. The contour levels for each contour image are automatically generated, but can be manually set with the `contour_levels0` keyword (see [**map_plot()**](https://github.com/kjdoore/spec_map_analysis/blob/master/spec_map_analysis/plots/map_plot.py) for more details). Finally, the color and opacity of each set of contours can be set with the `contour_color0` and `contour_alpha0` keywords, respectively. Both must be a list with `contour_color0` containing the colors of the respective sets of contours, and `contour_alpha0` containing the opacity scalars. 

For example, we can add white contours of the signal-to-noise to a line intensity image using:

```python
from spec_map_analysis.plots import map_plot

input_file = 'line_intensity.fits'
out_file = 'intensity_map_w_contours.png'
# Extension name of the image
extension = 'LINE_INTENSITY'

map_plot(out_file, input_file, extension, title0='Line Intensity Map', cbar_label0='W s$^{-2}$',
         contour0=[input_file], contour_extension0=['SNR'], contour_color0=['white'], contour_alpha0=[0.3])
```


### Multiple Images

Plotting multiple images with [**map_plot()**](https://github.com/kjdoore/spec_map_analysis/blob/master/spec_map_analysis/plots/map_plot.py) is relatively simple. Any additional images that are to be plotted in a new subplot can be read in using the generalized `input_fileX` and `extensionX` keywords, where the ending `X` is a proxy for a python style integer index indicating the image number.

For example if we had a FITS file containing the three image extensions of line intensity, line intensity uncertainty, and continuum emission for a fit emission line, we could plot all three images using:

```python
from spec_map_analysis.plots import map_plot

input_file = 'line_intensity.fits'
out_file = 'int_unc_cont_maps.png'

map_plot(out_file, input_file0=input_file, extension0='LINE_INTENSITY', 
         title0='Line Intensity', cbar_label0='W s$^{-2}$', 
         input_file1=input_file, extension1='UNCERTAINTY', 
         title1='Uncertainty', cbar_label1='W s$^{-2}$',
         input_file2=input_file, extension2='CONTINUUM', 
         title2='Continuum', cbar_label2='Jy')
```

This results in a saved graphic that has three subplots with each containing one of the images.

The shape of the subplot grid depends upon the number of input images. The grid will either be square or have more columns than rows. The exact shape is algorithmically generated, but the input images will always be plot in order from left to right, then top to bottom.

From the above example, it can be seen that a title and colorbar label can be given to each image. This is done using the generalized keywords `titleX` and `cbar_labelX`, where again the ending `X` is a proxy for a python style integer index indicating the image number. This generalization of keywords applies to all keywords ending in a `0` described above. To list them for convenience, it includes: `titleX`, `cbar_labelX`, `rangeX`, `apertureX` (and all of its properties), and `contourX` (and all of its properties).

**NOTE:** All plotted images will have the same axis range in terms of RA and Dec. The range of the 0th or first input image determines this axis range unless the `ralim` and `declim` keywords are set. In that case, all images will have the axis ranges provided in `ralim` and `declim`.

#### Image Outlines

When plotting multiple images that cover different sky areas, it may be convenient to add the outline of the smaller image onto the larger image. This can be done with the generalized `outlineX` keyword, which is set to the integer index of the image whose outline is to be over-plotted.

For example, we can plot two images that cover different sky areas, with the larger image being input first as to set the image axis range. We can then overlay the outline of the smaller image on both images to better see its location on the larger image using: 

```python
from spec_map_analysis.plots import map_plot

input_file0 = 'line_intensity_larger_area.fits'
input_file1 = 'line_intensity_smaller_area.fits'
out_file = 'line_intensity_maps.png'

map_plot(out_file, input_file0=input_file0, extension0='LINE_INTENSITY', 
         title0='Line Intensity', cbar_label0='W s$^{-2}$', outline0=1,
         input_file1=input_file1, extension1='LINE_INTENSITY', 
         title1='Line Intensity', cbar_label1='W s$^{-2}$', outline1=1)
```

Note that the `X` value in `outlineX` refers to the image that the outline will be place **ON**, and the value it is set to is the image whose outline is to be **DRAWN**.


### Performing Operations on Images

#### Operations with a Single Images

One of the most powerful components of [**map_plot()**](https://github.com/kjdoore/spec_map_analysis/blob/master/spec_map_analysis/plots/map_plot.py) is its ability to perform operations on an input image and plot the resulting image. This prevents the need to perform the operations outside of [**map_plot()**](https://github.com/kjdoore/spec_map_analysis/blob/master/spec_map_analysis/plots/map_plot.py), save the new image extension, and then input it into [**map_plot()**](https://github.com/kjdoore/spec_map_analysis/blob/master/spec_map_analysis/plots/map_plot.py). These operations are performed using the `operatorX` keyword, where again the ending `X` is a proxy for a python style integer index indicating the image number. The code evoked when this keyword is given makes use of the built-in python function, [eval()](https://docs.python.org/3/library/functions.html#eval). This function takes a string and performs the python statement within the string.

To give an example, we may want to plot the inverse of an image. This can be done using:

```python
from spec_map_analysis.plots import map_plot

input_file = 'data.fits'
out_file = 'inverse_map.png'

map_plot(out_file, input_file0=input_file, extension0=1, operator0=['1/a'])
```

Initially, it may be confusing why `operator0` is set to `'1/a'`. This is because the input image is read into [**map_plot()**](https://github.com/kjdoore/spec_map_analysis/blob/master/spec_map_analysis/plots/map_plot.py) as the variable `a` for convenience. Therefore, the [eval()](https://docs.python.org/3/library/functions.html#eval) function evaluates the statement `'1/a'` and returns the resulting inverted image.

If multiple operations are to be performed, they can be compressed into a single string, or each operation can be included as a separate string in the `operatorX` list. If multiple operations are given as separate strings, each operation is performed sequentially, with the output from each operation being saved as the variable `a`. For example, we can make two different calls to [**map_plot()**](https://github.com/kjdoore/spec_map_analysis/blob/master/spec_map_analysis/plots/map_plot.py), where we can take the square root of our squared input data (derive the absolute value), with each resulting in the same output image using:

```python
from spec_map_analysis.plots import map_plot

input_file = 'data.fits'
out_file = 'inverse_map.png'

# numpy is imported as np in map_plot
map_plot(out_file, input_file0=input_file, extension0=1, operator0=['np.sqrt(a**2)'])
# OR
map_plot(out_file, input_file0=input_file, extension0=1, operator0=['a**2', 'np.sqrt(a)'])
```

#### Operations with Multiple Images

Operations involving multiple images means that multiple images (all of which must have the same WCS information in the headers) will have operations performed on them to combine them into a single image for plotting. A common use for this would be to plot the ratio of data and uncertainty images (i.e., a signal-to-noise image), without having to make the image and save it to a separate image extension. Operations involving multiple images are similar to the single image case. However, extra images, besides the input image, will need to be specified using the keywords `operator_fileX` and `operator_extensionX`, where once again the ending `X` is a proxy for a python style integer index indicating the image number. Both keywords are lists with `operator_fileX` containing the name of the files to be read that will have an operation performed and `operator_extensionX` is the extension name or index. These images will be read in sequentially as the operations are performed as the variable `b`.

For example, to plot the SNR image as discussed above, we can use: 

```python
from spec_map_analysis.plots import map_plot

input_file = 'data.fits'
out_file = 'SNR_map.png'

map_plot(out_file, input_file0=input_file, extension0='DATA', operator0=['a/b'],
         operator_file0=[input_file], operator_extension0=['UNCERTAINTY'])
```

Here the secondary image (`operator_fileX`) is read in as the variable `b`, and divided by the original input, which is read in as the variable `a`.

If more than two images are to be combined with various operations, then multiple strings will need to be in the `operatorX` list. Each string will give the operation to be performed on the result of the previous operation, which is always saved as the variable `a`. For example, say we want to combine four different images as **w/x + y/z**, where each image is a unique letter. The `operatorX` keyword word would then need to be given as a list that preforms the required operations using only two images at a time. To do this, we would first need to rearrange the above expression to **(wz/x + y)/z**. Then we could do (assuming all image extensions are in the same file):

```python
from spec_map_analysis.plots import map_plot

input_file = 'data.fits'
out_file = 'complex_map.png'

map_plot(out_file, input_file0=input_file, extension0='W', operator0=['a*b', 'a/b', 'a+b', 'a/b'],
         operator_file0=[input_file, input_file, input_file, input_file], 
         operator_extension0=['Z', 'X', 'Y', 'Z'], title0='w/x+y/x=(wz/x+y)/z')
```

Note that the `operator_fileX` first image extension is read in as `b`, then the first operation in the `operatorX` list is performed and the result save to the variable `a`. Next the second image is read in as `b` as well, then the second operation in the `operatorX` list is performed, and the result again is saved to the variable `a`. This is repeated until all operations are performed and shows how complex operations with several images are performed. 

Generalized, an image in the `operator_fileX` and `operator_extensionX` list is read in as `b`, the corresponding indexed operation string in `operatorX` is performed, and the operation output is save as `a`. Highly complex operations are allowed when using `operatorX`. However, each operation string in the list can only contain at most two images, not reference any other variables besides `a` and `b`, and only use built-in python functions and numpy functions (in the standard `import numpy as np` fashion).

**NOTE:** If an operation is performed only on `a`, it does not require a corresponding index in the `operator_fileX` and `operator_extensionX` lists if it is the last operation or only operation performed.