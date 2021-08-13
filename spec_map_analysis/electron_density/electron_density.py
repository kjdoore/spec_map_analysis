def electron_density(out_filename, filename_line_ratio, element, spec, temperature, aperture=False,
                     aperture_only=False, **kwargs):
    """
    A function that reads in a FITS file with a line ratio
      map. The map is converted into an electron density map
      using PyNeb and the specified temperature.
      If aperture photometry was performed and line ratios
      for the apertures are given in a table, then electron
      densities and determined for each aperture and are
      appended to the table extension on the input FITS file,
      and included in the output FITS file.
    Parameters
    ----------
    out_filename : string
         A string containing the name of the resulting electron
         density FITS file.
    filename_line_ratio : string
        A string containing the line ratio FITS file to be read.
    element : string
        A string containing the atomic symbol for the line ratio
        element.
    spec : integer
        An integer containing the ionization stage of the line ratio
        element.
    temperature : scalar
        A scalar containing the temperature at which to determine
        the electron density of the line ratio in Kelvin.
    aperture : boolean, optional
        A boolean stating whether aperture photometry and line
        ratios are given in a table extension, and electron
        densities should be computed for each aperture.
    aperture_only : boolean, optional
        A boolean stating if electron densities should only
        be calculated for the aperture photometry and
        line ratios given in a table extension, and not both
        the apertures and map.
    kwargs
        Keyword arguments passed to PyNeb.getTemDen
    """

    import pyneb as pn
    from astropy.io import fits
    from astropy.table import Table
    import numpy as np

    # Read in the line ratio file
    file = fits.open(filename_line_ratio)

    # Generate the ion given the element and species
    ion = pn.Atom(element, spec)

    # Do if electron densities maps are wanted
    if not aperture_only:
        # Place the data and uncertainty into arrays
        int_ratio = file['LINE_RATIO'].data
        int_ratio_unc = file['UNCERTAINTY'].data

        # Compute the electron densities for the line ratio map
        density = ion.getTemDen(int_ratio, tem=temperature, **kwargs)
        # Compute the upper and lower uncertainties for the density
        density_plus = ion.getTemDen(int_ratio + int_ratio_unc, tem=temperature, **kwargs)
        density_minus = ion.getTemDen(int_ratio - int_ratio_unc, tem=temperature, **kwargs)

        # Determine if density_plus or density_minus is the upper or lower uncertainty
        #   (i.e., larger line ratios imply lower densities if the lower energy
        #   transition is given as the numerator in the line ratio)
        if np.nansum(density_plus - density) > 0:
            density_upper = density_plus - density
            density_lower = density - density_minus
        else:
            density_upper = density_minus - density
            density_lower = density - density_plus

    # Compute the electron densities of the aperture line ratios and add them to the FITS table
    if aperture:
        aper_tbl = Table(file['APERTURE_LINE_RATIO'].data)
        aper_int_ratio = aper_tbl['line_ratio']
        aper_int_ratio_unc = aper_tbl['line_ratio_unc']
        aper_density = ion.getTemDen(aper_int_ratio, tem=temperature, **kwargs)
        aper_tbl.add_column(aper_density, name='aperture_electron_density')
        aper_density_plus = ion.getTemDen(aper_int_ratio + aper_int_ratio_unc, tem=temperature, **kwargs)
        aper_density_minus = ion.getTemDen(aper_int_ratio - aper_int_ratio_unc, tem=temperature, **kwargs)
        if np.nansum(aper_density_plus - aper_density) > 0:
            aper_tbl.add_column(aper_density_plus - aper_density, name='aperture_electron_density_upper_unc')
            aper_tbl.add_column(aper_density - aper_density_minus, name='aperture_electron_density_lower_unc')
        else:
            aper_tbl.add_column(aper_density_minus - aper_density, name='aperture_electron_density_upper_unc')
            aper_tbl.add_column(aper_density - aper_density_plus, name='aperture_electron_density_lower_unc')

    # Generate the image extensions with the input extension header updated
    if not aperture_only:
        image_hdr = file[1].header
        image_hdr['BUNIT'] = 'n_e/pixel'
        image_hdr['EXTNAME'] = 'ELECTRON_DENSITY'
        el_den_hdu = fits.ImageHDU(density, header=image_hdr)

        image_hdr = file[1].header
        image_hdr['BUNIT'] = 'n_e/pixel'
        image_hdr['EXTNAME'] = 'LOWER_UNCERTAINTY'
        el_den_unc_hdu_low = fits.ImageHDU(density_lower, header=image_hdr)

        image_hdr = file[1].header
        image_hdr['BUNIT'] = 'n_e/pixel'
        image_hdr['EXTNAME'] = 'UPPER_UNCERTAINTY'
        el_den_unc_hdu_up = fits.ImageHDU(density_upper, header=image_hdr)

    # Add a HISTORY to the primary header stating the parameters used in the
    #   electron density calculation
    primary_hdr = file[0].header
    primary_hdr['HISTORY'] = 'Electron Density Calculation'
    primary_hdr['HISTORY'] = '-- Electron density was calculated using PyNeb'
    primary_hdr['HISTORY'] = '-- Parameter inputs for electron density Calculation'
    primary_hdr['HISTORY'] = '  element = ' + element
    primary_hdr['HISTORY'] = '  species = ' + str(spec)
    primary_hdr['HISTORY'] = '  temperature = ' + str(temperature) + ' K'
    for key in kwargs:
        primary_hdr['HISTORY'] = '  ' + key + ' = ' + str(kwargs[key])
    primary_hdr['HISTORY'] = '--'
    primary_hdr['HISTORY'] = ''

    primary_hdu = fits.PrimaryHDU(header=primary_hdr)

    # If aperture electron densities were calculated, place them in a FITS table, and
    #   generate FITS HDU.
    if aperture:
        table_hdu = fits.BinTableHDU(data=aper_tbl, name='APERTURE_ELECTRON_DENSITY')
        if aperture_only:
            hdul_new = fits.HDUList([primary_hdu, table_hdu])
        else:
            hdul_new = fits.HDUList([primary_hdu, el_den_hdu, el_den_unc_hdu_low, el_den_unc_hdu_up, table_hdu])
    else:
        hdul_new = fits.HDUList([primary_hdu, el_den_hdu, el_den_unc_hdu_low, el_den_unc_hdu_up])

    hdul_new.writeto(out_filename, overwrite=True)
