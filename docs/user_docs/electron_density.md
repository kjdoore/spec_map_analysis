# Line Ratios and Electron Densities ([spec_map_analysis.electron_density](https://github.com/kjdoore/spec_map_analysis/blob/master/docs/user_docs/electron_density.md))

## Introduction

In **spec_map_analysis**, the **electron_density** sub-package is used to derive line ratio maps of two line intensity images and electron density maps from the line ratio maps. The electron density maps are derived from line ratios using the [PyNeb](http://research.iac.es/proyecto/PyNeb//) atomic library. Aperture values for both line ratios and electron densities can optionally be computed from aperture photometry of the two input line intensity images. The base result is a FITS file for both the line ratio maps and electron density maps with extensions containing the data and uncertainty maps.


## Line Ratio Maps

To only generate a line ratio map and its corresponding uncertainty map, a simple call to the [**line_ratio()**](https://github.com/kjdoore/spec_map_analysis/blob/master/spec_map_analysis/electron_density/line_ratio.py) function is required. The required inputs are the name of the output file, the names of the two line intensity input files, and the intensity files' extensions for the data and uncertainty images. The function simply takes the line intensity and uncertainty maps, which must be convolved and resampled to the same pixel scale and spatial resolution (preferable using [**spec_map_analysis.reproject_convolve.reproject_convolve()**](https://github.com/kjdoore/spec_map_analysis/blob/master/docs/user_docs/reproject_convolve.md)), and determines the ratio and propagates uncertainties. The first input file will be considered the numerator in the line ratio, and the second file will be the denominator.

For an example, we can derive the line ratio and uncertainties using the convolved [OIII] 52&mu;m line intensity maps and resampled [OIII] 88&mu;m line intensity maps from the [**spec_map_analysis.reproject_convolve.reproject_convolve()**](https://github.com/kjdoore/spec_map_analysis/blob/master/docs/user_docs/reproject_convolve.md) example using:

```python
from spec_map_analysis.electron_density import line_ratio

# Line intensity will be the numerator in the ratio 
input_file1 = 'line_intensity_03_88_resampled.fits'
# Extension names of the 1st input file
extensions1 = ['LINE_INTENSITY', 'UNCERTAINTY']

# Line intensity will be the denominator in the ratio 
input_file2 = 'line_intensity_03_52_convolved.fits'
# Extension names of the 2nd input file
extensions2 = ['LINE_INTENSITY', 'UNCERTAINTY']

output_file = 'line_ratio_03_88_52.fits'

line_ratio(output_file, input_file1, input_file2, extensions1, extensions2)
```

This call to [**line_ratio()**](https://github.com/kjdoore/spec_map_analysis/blob/master/spec_map_analysis/electron_density/line_ratio.py) saves the line ratio map and uncertainty map as image extensions in the output file.


## Electron Density Maps

Electron density maps can be derived using the [**electron_density()**](https://github.com/kjdoore/spec_map_analysis/blob/master/spec_map_analysis/electron_density/electron_density.py) function. The function requires five arguments: the name of the output file, the name of the input line ratio file, the element of the emission lines, the ionization stage of said element, and the temperature of the electrons in **K**. Also, one or more keywords need to be given specifying which lines are used in the input line ratio (see [getTemDen()](https://morisset.github.io/PyNeb_Manual/html/classpyneb_1_1core_1_1pynebcore_1_1_atom.html#ad245cc7a9a33c7db1c2276b53ca02537) for details).

For example, we can determine the electron density of the line ratio in the above example assuming an electron temperature of 10<sup>4</sup> K using:

```python
from spec_map_analysis.electron_density import electron_density

# Line ratio file from line_ratio()
input_file = 'line_ratio_03_88_52.fits'
output_file = 'electron_density.fits'

# Parameters for PyNeb
element = 'O'
species = 3
temperature = 1e4

electron_density(output_file, input_file, element, species, temperature, 
                 to_eval='L(88e4)/ L(52e4)')
```

Note the `to_eval` keyword specifies which [OIII] line ratio is used. We could have also equally have used `wave1` and `wave2`; or `lev_j1`, `lev_j2`, `lev_i1`, and `lev_i2` as detailed in [getTemDen()](https://morisset.github.io/PyNeb_Manual/html/classpyneb_1_1core_1_1pynebcore_1_1_atom.html#ad245cc7a9a33c7db1c2276b53ca02537). 

The resulting output file contains three image extensions. One is the electron density map, and the others are the upper and lower uncertainty ranges since the conversion from line ratio to electron density is not a symmetric relation.



## Aperture Line Ratios and Electron Densities

### Line Ratios
In addition to maps of line ratios and electron densities, aperture measurements can also be calculated on the line intensity maps while computing line ratios in the [**line_ratio()**](https://github.com/kjdoore/spec_map_analysis/blob/master/spec_map_analysis/electron_density/line_ratio.py) function. To compute aperture measurements, the keywords `aper_pos` and `aper_size` need to be specified. `aper_pos` is a list of 2 element tuples specifying the position of any aperture measurements to be taken on the image. The units of the positions can be specified in pixel or sky units depending on the value of of the `aper_unit` keyword (`'sky'` (default) or `'pixel'`). If the position are specified in sky units, they must be a form acceptable by [astropy.coordinates.SkyCoord](https://docs.astropy.org/en/stable/api/astropy.coordinates.SkyCoord.html#astropy.coordinates.SkyCoord) in terms of Right Ascension and Declination. `aper_size` is a list of tuples giving the size of the aperture in **arcseconds** if `aper_unit='sky'` or in **pixels** if `aper_unit='pixel'` and the tuple length depends on the aperture shape. The shape of the aperture is specified with `aper_shape`, which can be either `'circular'` (default), `'elliptical'`, or `'rectangular'`.

For example, we can derive circular aperture measurements at three different locations in the above line ratio example. Here we will assume the example is of NGC 1569, and we will take two apertures of different sizes at the central peak and another at the secondary emission peak. To generate line ratio values at these apertures with radii of  6'' and 12'', we can use:

```python
from spec_map_analysis.electron_density import line_ratio

# Line intensity will be the numerator in the ratio 
input_file1 = 'line_intensity_03_88_resampled.fits'
# Extension names of the 1st input file
extensions1 = ['LINE_INTENSITY', 'UNCERTAINTY']

# Line intensity will be the denominator in the ratio 
input_file2 = 'line_intensity_03_52_convolved.fits'
# Extension names of the 2nd input file
extensions2 = ['LINE_INTENSITY', 'UNCERTAINTY']

output_file = 'line_ratio_03_88_52_with_aperture.fits'

line_ratio(output_file, input_file1, input_file2, extensions1, extensions2, 
           aper_pos=[('4h30m47.0456s', '64d51m00.506s'), ('4h30m47.0456s', '64d51m00.506s'), 
                     ('4h30m52.1933s', '64d50m46.827s')],  aper_size=[(12) ,(6), (6)])
```

The result will contain, in addition to the line ratio and uncertainty maps, a table containing aperture measurements and the corresponding line ratios.

If only the table with the aperture measurements is wanted, then the `aperture_only` keyword can be set to `True`, which suppresses the calculation of the line ratio map.

**Note**: See the description of [**line_ratio()**](https://github.com/kjdoore/spec_map_analysis/blob/master/spec_map_analysis/electron_density/line_ratio.py) for a full list and detailed description of all other parameters related to apertures not discussed here.

### Electron Densities

Computing aperture electron densities from the aperture line ratios is simple. Using the [**electron_density()**](https://github.com/kjdoore/spec_map_analysis/blob/master/spec_map_analysis/electron_density/electron_density.py) function, if the input line ratio file contains aperture line ratios derived by the [**line_ratio()**](https://github.com/kjdoore/spec_map_analysis/blob/master/spec_map_analysis/electron_density/line_ratio.py) function, simply set the keyword `aperture=True`. This tells the function to expect a table with aperture line ratio values and to compute the corresponding electron densities. Just like the [**line_ratio()**](https://github.com/kjdoore/spec_map_analysis/blob/master/spec_map_analysis/electron_density/line_ratio.py) function, if only the table with the aperture measurements is wanted, the `aperture_only` keyword can be set to `True`, resulting in a file with only an aperture table extension.
