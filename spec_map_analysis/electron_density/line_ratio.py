def line_ratio(out_filename, filename1, filename2, extensions1, extensions2,
               aper_pos=None, aper_unit='sky', aper_shape='circular', aper_size=None, aper_mask=None,
               aper_snr_mask=None, aperture_only=False, mc_iter=100, **kwargs):
    """
    A function that reads two FITS files with line intensity
      maps. The two maps, which need to be aligned and on the
      same pixel scale, are combined into a line ratio map,
      and saved to an output FITS file.
      There is the option to preform aperture photometry on
      the line intensity maps using input positions and
      aperture shapes. This photometry is then used to
      compute line ratios and placed in a table extension
      in the output FITS file.
    Parameters
    ----------
    out_filename : string
         A string containing the name of the resulting line
         ratio FITS file.
    filename1 : string
        A string containing the first line intensity FITS file
        to be read. This file is the numerator of the line ratio.
    filename2 : string
        A string containing the first line intensity FITS file
        to be read. This file is the denominator of the line ratio.
    extensions1 : two element list of strings or integers
        A two element list of string and/or integers containing the name
        or index of the extensions to be read for filename1. The order
        of the list must be 0) the data and 1) the uncertainty of the data.
    extensions2 : two element list of strings or integers
        A two element list of string and/or integers containing the name
        or index of the extensions to be read for filename2. The order
        of the list must be 0) the data and 1) the uncertainty of the data.
    aper_pos : list of tuples, optional
        A list of 2 element tuples containing the position(s) of
        the aperture(s) in units specified by aper_unit. If
        aper_unit='sky', then position(s) are specified in sky units,
        and must be a form acceptable by astropy.coordinates.SkyCoord.
        The default reference frame is ICRS, but it can be changed
        with the 'frame' additional keyword argument. If
        aper_unit='pixel', then position(s) are specified in pixel units.
    aper_unit : string, optional
        A case-insensitive string containing the unit for the aperture
        position and size. The options are 'sky' (default) and 'pixel'.
        'sky' units refer to sky units as in an angle in reference to
        a given sky coordinate system, and 'pixel' units refer to the
        native pixel units of the map array. The unit will be the same
        for all apertures.
    aper_shape : string, optional
        A case-insensitive string specifying the shape of the aperture.
        The options are 'circular' (default), 'elliptical', and
        'rectangular'. The aperture shape will be the same for all apertures.
    aper_size : list of tuples, optional
        A list of 1 (3) element tuples containing the size (and orientation)
        of the circular (elliptical or rectangular) apertures. For circular
        aperture(s), the 1 element tuple(s) contain the radius of the
        aperture. For elliptical aperture(s), the 3 element tuple(s)
        contain the semimajor, semiminor, and position angle of the
        ellipse(s), respectively. For rectangular aperture(s), the 3
        element tuple(s) contain the full width, full height, and position
        angle of the rectangle(s), respectively. If aper_unit='sky', all
        sizes must be given in arcseconds, and if aper_unit='pixel', all
        sizes must be given in pixels. Regardless of aper_unit, position
        angles must be given in radians.
    aper_mask : 2D boolean array-like, optional
        A two dimensional boolean array with the same shape as the input
        file data. A True value indicates the corresponding element
        of the data that are to excluded from all aperture calculations.
    aper_snr_mask : scalar, optional
        A scalar indicating if a mask for the aperture calculations
        should be made that limits the line intensities of the
        numerator and denominator to pixels of the line ratio
        that have a SNR greater than or equal to the given value.
        Pixels with a line ratio SNR less than the value are
        excluded from all aperture calculations. Note: Cannot
        be set if aperture_only=True.
    aperture_only : boolean, optional
        A boolean stating if only aperture line ratios should
        be calculated from the input line intensity maps, and
        not both the aperture line ratios and a line ratio map.
    mc_iter : integer, optional
        An integer giving the number of Monte Carlo (MC) iterations
        to run in computing the aperture photometry and uncertainty.
    kwargs
        Keyword arguments passed to astropy.coordinates.SkyCoord
    """

    from astropy.io import fits
    import numpy as np
    from photutils import aperture
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    from astropy.wcs import WCS
    from astropy.table import Table

    # Read in the line intensity files
    numerator_file = fits.open(filename1)
    denominator_file = fits.open(filename2)


    if not aperture_only:
        # Determine the line ratio
        line_ratio = numerator_file[extensions1[0]].data/denominator_file[extensions2[0]].data
        # Propagate the line intensity uncertainty to get the line ratio uncertainty, assuming no covariance
        line_ratio_unc = line_ratio * np.sqrt((numerator_file[extensions1[1]].data/
                                               numerator_file[extensions1[0]].data)**2 +
                                              (denominator_file[extensions2[1]].data/
                                               denominator_file[extensions2[0]].data)**2)

    # Determine the total line intensity within each aperture if specified.
    if aper_pos is not None:
        # Create and combine masks needed for aperture calculations
        # Create mask for NaN values
        mask = np.logical_or(np.isnan(numerator_file[extensions1[0]].data),
                             np.isnan(denominator_file[extensions2[0]].data))
        # Create mask for line ratio SNR values if given, and combine with NaN mask
        if aper_snr_mask is not None:
            line_ratio_snr = line_ratio/line_ratio_unc
            snr_mask = line_ratio_snr < aper_snr_mask
            mask = np.logical_or(snr_mask, mask)
        # Combine NaN/SNR mask with input mask if given
        if aper_mask is not None:
            mask = np.logical_or(aper_mask, mask)

        aper_numerator = np.zeros(len(aper_pos))
        aper_numerator_unc = np.zeros(len(aper_pos))
        aper_denominator = np.zeros(len(aper_pos))
        aper_denominator_unc = np.zeros(len(aper_pos))

        for i in range(0, len(aper_pos)):
            # Use the specified aperture shape and pixel or sky unit
            if aper_shape.casefold() == 'circular':
                if aper_unit.casefold() == 'sky':
                    position = SkyCoord([aper_pos[i]], **kwargs)
                    aper = aperture.SkyCircularAperture(position, aper_size[i] * u.arcsec)
                elif aper_unit.casefold() == 'pixel':
                    aper = aperture.CircularAperture(aper_pos[i], aper_size[i])
                else:
                    raise ValueError("aper_unit must be 'sky' or 'pixel'")
            elif aper_shape.casefold() == 'elliptical':
                if aper_unit.casefold() == 'sky':
                    position = SkyCoord([aper_pos[i]], **kwargs)
                    aper = aperture.SkyEllipticalAperture(position, aper_size[i][0] * u.arcsec,
                                                          aper_size[i][1] * u.arcsec,
                                                          theta=aper_size[i][2] * u.rad)
                elif aper_unit.casefold() == 'pixel':
                    aper = aperture.EllipticalAperture(aper_pos[i], aper_size[i][0], aper_size[i][1],
                                                       theta=aper_size[i][2])
                else:
                    raise ValueError("aper_unit must be 'sky' or 'pixel'")
            elif aper_shape.casefold() == 'rectangular':
                if aper_unit.casefold() == 'sky':
                    position = SkyCoord([aper_pos[i]], **kwargs)
                    aper = aperture.SkyRectangularAperture(position, aper_size[i][0] * u.arcsec,
                                                           aper_size[i][1] * u.arcsec,
                                                           theta=aper_size[i][2] * u.rad)
                elif aper_unit.casefold() == 'pixel':
                    aper = aperture.RectangularAperture(aper_pos[i], aper_size[i][0], aper_size[i][1],
                                                        theta=aper_size[i][2])
                else:
                    raise ValueError("aper_unit must be 'sky' or 'pixel'")
            else:
                raise ValueError("aper_shape must be 'circular', 'elliptical', or 'rectangular'")

            # Compute the aperture photometry and uncertainties depending on the aperture units
            #   for both line maps, using a Monte Carlo method.
            mc_numer = np.zeros(mc_iter)
            mc_denom = np.zeros(mc_iter)
            for j in range(0, mc_iter):
                # perturb data by the normally distributed uncertainty
                numer_data = numerator_file[extensions1[0]].data + numerator_file[extensions1[1]].data * \
                             np.random.standard_normal(numerator_file[extensions1[1]].data.shape)
                denom_data = denominator_file[extensions2[0]].data + denominator_file[extensions2[1]].data * \
                             np.random.standard_normal(denominator_file[extensions2[1]].data.shape)
                if aper_unit.casefold() == 'sky':
                    numer_aper_phot = aperture.aperture_photometry(numer_data, aper, mask=mask,
                                                                   wcs=WCS(numerator_file[extensions1[0]].header))
                    denom_aper_phot = aperture.aperture_photometry(denom_data, aper, mask=mask,
                                                                   wcs=WCS(denominator_file[extensions2[0]].header))
                elif aper_unit.casefold() == 'pixel':
                    numer_aper_phot = aperture.aperture_photometry(numer_data, aper, mask=mask)
                    denom_aper_phot = aperture.aperture_photometry(denom_data, aper, mask=mask)
                else:
                    raise ValueError("aper_unit must be 'sky' or 'pixel'")

                mc_numer[j] = numer_aper_phot[0]['aperture_sum']
                mc_denom[j] = denom_aper_phot[0]['aperture_sum']

            aper_numerator[i] = np.median(mc_numer)
            aper_numerator_unc[i] = np.std(mc_numer)
            aper_denominator[i] = np.median(mc_denom)
            aper_denominator_unc[i] = np.std(mc_denom)

        # Determine the line ratios for each aperture
        aper_line_ratio = aper_numerator / aper_denominator
        aper_line_ratio_unc = aper_line_ratio * np.sqrt((aper_numerator_unc / aper_numerator)**2 +
                                                        (aper_denominator_unc / aper_denominator)**2)

        # Place results and parameters into a table with units
        table_arr = {'aperture_size': aper_size, 'aperture_center': aper_pos,
                     'aperture_unit': [aper_unit] * len(aper_pos), 'aperture_shape': [aper_shape] * len(aper_pos),
                     'numerator_aperture_sum': aper_numerator, 'numerator_aperture_sum_unc': aper_numerator_unc,
                     'denominator_aperture_sum': aper_denominator, 'denominator_aperture_sum_unc': aper_denominator_unc,
                     'line_ratio': aper_line_ratio, 'line_ratio_unc': aper_line_ratio_unc}
        if aper_unit.casefold() == 'sky':
            unit = ['arcsec', None, None, None, 'W/m^2', 'W/m^2', 'W/m^2', 'W/m^2', None, None]
        elif aper_unit.casefold() == 'pixel':
            unit = ['pixel', 'pixel', None, None, 'W/m^2', 'W/m^2', 'W/m^2', 'W/m^2', None, None]
        else:
            raise ValueError("aper_unit must be 'sky' or 'pixel'")

        tablehdu = Table(data=table_arr, units=unit)

    if not aperture_only:
        # Change the image extension header
        image_hdr = numerator_file[extensions1[0]].header
        image_hdr['BUNIT'] = 'ratio/pixel'
        image_hdr['EXTNAME'] = 'LINE_RATIO'
        line_ratio_hdu = fits.ImageHDU(line_ratio, header=image_hdr)

        image_hdr = numerator_file[extensions1[1]].header
        image_hdr['BUNIT'] = 'ratio/pixel'
        image_hdr['EXTNAME'] = 'UNCERTAINTY'
        line_ratio_unc_hdu = fits.ImageHDU(line_ratio_unc, header=image_hdr)

    # Generate the primary header
    primary_hdr = numerator_file[0].header
    # Add a HISTORY to the primary header if aperture masks were used
    if aper_mask is not None or aper_snr_mask is not None:
        primary_hdr['HISTORY'] = 'Masks Used in Aperture Photometry'
        primary_hdr['HISTORY'] = '-- Aperture Mask Parameters'
        if aper_snr_mask is not None:
            primary_hdr['HISTORY'] = '  aper_mask = True'
        if aper_snr_mask is not None:
            primary_hdr['HISTORY'] = '  aper_snr_mask = ' + str(aper_snr_mask)
        primary_hdr['HISTORY'] = '--'
        primary_hdr['HISTORY'] = ''
    primary1_hdu = fits.PrimaryHDU(header=primary_hdr)

    # If apertures were used, create table extension
    if aper_pos is not None:
        table_hdu = fits.BinTableHDU(data=tablehdu, name='APERTURE_LINE_RATIO')

    # Create the final HDU list
    if aperture_only:
        hdul_new = fits.HDUList([primary1_hdu, table_hdu])
    elif aper_pos is not None:
        hdul_new = fits.HDUList([primary1_hdu, line_ratio_hdu, line_ratio_unc_hdu, table_hdu])
    else:
        hdul_new = fits.HDUList([primary1_hdu, line_ratio_hdu, line_ratio_unc_hdu])
    hdul_new.writeto(out_filename, overwrite=True)
