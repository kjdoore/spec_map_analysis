def reproject_convolve(out_filename1, out_filename2, filename1, filename2, spatial_res1, spatial_res2,
                       file_to_reproject=2, reprojection_method='interp', data_extension='LINE_INTENSITY',
                       unc_extension='UNCERTAINTY', fwhm=False, pixel_unit='pixel', **kwargs):
    """
    A function that reads two input FITS files and reprojects (resamples)
      the the specified image, and then convolves the first image with a
      Gaussian of the specified sigma (sigma_smooth). The uncertainties
      of the convolved image are propagated assuming the uncertainties are
      correlated like a Gaussian convolution with the given beam size.
      The resulting images are then saved as FITS image extensions to
      the specified output files.
    Parameters
    ----------
    out_filename1 : string
        A string containing the file name for filename1 to be writen
        to after being reprojected and/or convolved
    out_filename2 : string
        A string containing the file name for filename2 to be writen
        to after being reprojected
    filename1 : string
        A string containing the file to be opened, which has a higher
        spatial resolution than filename2 and will be convolved with
        a Gaussian of the given sigma, sigma_smooth
    filename2 : string
        A string containing the secondary file to be opened, which
        has a lower spatial resolution than filename1 and
        can be reprojected (resampled) to the filename1 pixel scale
    spatial_res1 : scalar
        A scalar giving the spatial resolution of the first image,
        which will be convolved, in arcseconds. The value must be
        the standard deviation of the spatial resolution, unless
        the keyword fwhm is set to True.
    spatial_res2 : scalar
        A scalar giving the spatial resolution of the second image
        in arcseconds. The value must be the standard deviation of
        the spatial resolution, unless the keyword fwhm is set to True.
    file_to_reproject : integer, optional
        A value of 1 or 2, which indicates which input file is to be
        reprojected (resampled)
    reprojection_method : string, optional
        A case-insensitive string containing the reprojection method to use.
        Can be any of the following: 'adaptive', 'exact', or 'interp' (default).
        See https://reproject.readthedocs.io/en/stable/celestial.html for
        details on the reproject functions.
    data_extension : string or int, optional
        A string or integer containing the name or index of the data extension
        for both input files (both files must have the same extension name or
        index)
    unc_extension : string or int, optional
        A string or integer containing the name or index of the uncertainty
        extension for both input files (both files must have the same
        extension name or index)
    fwhm : boolean
        A boolean stating if the spatial resolutions are given as FWHMs (if
        True) or standard deviations (if False)
    pixel_unit : string, optional
        A case insensitive string containing the type of unit in each
        pixel. For example, a flux in Jy for each pixel could be expressed
        in Jy/pixel (the pixel unit would be per pixel) or Jy/sr (the pixel
        unit would be per anglular unit). This is needed to properly scale the
        reprojected (resampled) image if the pixel units are per pixel.
        Available inputs are: 'pixel' (default) and 'angular'
    kwargs
        Keyword arguments passed to astropy.convolution.GaussianXDKernel
    """

    import numpy as np
    from astropy.io import fits
    from spec_map_analysis.reproject_convolve import gauss_convol_with_corr_error as gcce
    import reproject as rpjct

    # Read in the files depending on which image is to be reprojected (resampled)
    if file_to_reproject == 1:
        hdu_reprj = fits.open(filename1)
        hdu_noreprj = fits.open(filename2)
    elif file_to_reproject == 2:
        hdu_reprj = fits.open(filename2)
        hdu_noreprj = fits.open(filename1)
    else:
        raise ValueError("file_to_reproject must have a value of 1 or 2")

    # Convert the input spatial resolutions of the input images into sigma (if necessary)
    #   and compute the sigma needed for the Gaussian kernel.
    if fwhm:
        sigma_beam = spatial_res1 / (2 * np.sqrt(2 * np.log(2)))
        sigma_smooth = np.sqrt(spatial_res2 ** 2 - spatial_res1 ** 2) / (2 * np.sqrt(2 * np.log(2)))
    else:
        sigma_beam = spatial_res1
        sigma_smooth = np.sqrt(spatial_res2 ** 2 - spatial_res1 ** 2)

    # Reproject (resample) the specified image using the given method. Both the data and uncertainty
    #   need to be resampled.
    if reprojection_method.casefold() == 'adaptive':
        data_reprj = rpjct.reproject_adaptive(hdu_reprj, hdu_noreprj[data_extension].header,
                                              hdu_in=data_extension, return_footprint=False)
        unc_reprj = rpjct.reproject_adaptive(hdu_reprj, hdu_noreprj[unc_extension].header,
                                             hdu_in=unc_extension, return_footprint=False)
    elif reprojection_method.casefold() == 'exact':
        data_reprj = rpjct.reproject_exact(hdu_reprj, hdu_noreprj[data_extension].header,
                                           hdu_in=data_extension, return_footprint=False)
        unc_reprj = rpjct.reproject_exact(hdu_reprj, hdu_noreprj[unc_extension].header,
                                          hdu_in=unc_extension, return_footprint=False)
    elif reprojection_method.casefold() == 'interp':
        data_reprj = rpjct.reproject_interp(hdu_reprj, hdu_noreprj[data_extension].header,
                                            hdu_in=data_extension, return_footprint=False)
        unc_reprj = rpjct.reproject_interp(hdu_reprj, hdu_noreprj[unc_extension].header,
                                           hdu_in=unc_extension, return_footprint=False)
    else:
        raise ValueError("reprojection_method must be either 'interp'," +
                         "'exact', or 'adaptive'")

    # Rescale the reprojected (resampled) images based on the pixel_unit input
    # Assuming CDELT1 is in degrees as is standard FITS format
    reprj_hd = hdu_reprj[data_extension].header
    noreprj_hd = hdu_noreprj[data_extension].header
    if pixel_unit.casefold() == 'pixel':
        # Divide by scale of original sampling, multiply by new scale
        scale_factor = np.abs((noreprj_hd['CDELT1']*noreprj_hd['CDELT2'])/(reprj_hd['CDELT1']*reprj_hd['CDELT2']))
    elif pixel_unit.casefold() == 'angular':
        scale_factor = 1
    else:
        raise ValueError("pixel_area_unit_type must be either 'pixel' or 'angular'")
    data_reprj = data_reprj * scale_factor
    unc_reprj = unc_reprj * scale_factor

    # Prepare the first input image data and uncertainty to be convolved
    if file_to_reproject == 1:
        # If the first image was resampled, then convolve the resampled data
        data = data_reprj
        unc = unc_reprj

        # Convert sigma_smooth and sigma_beam to pixels
        # CDELTn in degrees per pixel
        # Check to make sure the pixels scale is the same for both data and uncertainty images
        # Use CDELTn from the non-resampled image, as the data currently has that pixel scale
        hd_data = hdu_noreprj[data_extension].header
        hd_unc = hdu_noreprj[unc_extension].header
        cdelt_data = hd_data['CDELT*']
        cdelt_unc = hd_unc['CDELT*']
        if np.abs(cdelt_data[0]) == np.abs(cdelt_unc[0]):
            sigma_smooth_pix = sigma_smooth / 3.6e3 / np.abs(cdelt_data[0])
            sigma_beam_pix = sigma_beam / 3.6e3 / np.abs(cdelt_data[0])
        else:
            raise ValueError('pixel scale (CDELTn) between images is not the same')
    else:
        # If the second image was resampled, then convolve the first image which was not resampled
        data = hdu_noreprj[data_extension].data
        unc = hdu_noreprj[unc_extension].data

        # Convert sigma_smooth and sigma_beam to pixels
        # CDELTn in degrees per pixel
        # Check to make sure the pixels scale is the same for both data and uncertainty images
        # Use CDELTn from the non-resampled image, as the first image currently has that pixel scale
        hd_data = hdu_noreprj[data_extension].header
        hd_unc = hdu_noreprj[unc_extension].header
        cdelt_data = hd_data['CDELT*']
        cdelt_unc = hd_unc['CDELT*']
        if np.abs(cdelt_data[0]) == np.abs(cdelt_unc[0]):
            sigma_smooth_pix = sigma_smooth / 3.6e3 / np.abs(cdelt_data[0])
            sigma_beam_pix = sigma_beam / 3.6e3 / np.abs(cdelt_data[0])
        else:
            raise ValueError('pixel scale (CDELTn) between images is not the same')

    # Convolve the first input image
    convol_data, convol_unc = gcce(data, unc, sigma_smooth_pix, sigma_beam_pix, **kwargs)

    # Create the headers for final data products. Edit the Primary headers to note the resampling
    #   and convolution in HISTORY comments. Also, update the image extension headers to reflect this,
    #   but only include the HISTORY in the primary
    if file_to_reproject == 1:
        # If the first file is reprojected (resampled) then state it was both reprojected
        #   and convolved in the primary header HISTORY
        primary_hdr1 = hdu_reprj[0].header
        primary_hdr1['HISTORY'] = 'Reprojection and Convolution'
        primary_hdr1['HISTORY'] = '-- Map was reprojected (resampled) and convolved'
        primary_hdr1['HISTORY'] = '-- Parameter inputs'
        primary_hdr1['HISTORY'] = '  sigma_smooth = ' + str(sigma_smooth) + ' arcsec'
        primary_hdr1['HISTORY'] = '  sigma_beam = ' + str(sigma_beam) + ' arcsec'
        primary_hdr1['HISTORY'] = '  reprojection_method = ' + str(reprojection_method)
        primary_hdr1['HISTORY'] = '  pixel_area_unit_type = ' + pixel_unit
        for key in kwargs:
            primary_hdr1['HISTORY'] = '  ' + key + ' = ' + str(kwargs[key])
        primary_hdr1['HISTORY'] = '--'
        primary_hdr1['HISTORY'] = ''
        # If the header had an extension name, remove it
        if 'EXTNAME' in primary_hdr1:
            primary_hdr1.remove('EXTNAME')
        # If a keyword is in both the primary and non-reprojected image header, replace
        #   the primary value with the non-reprojected value
        for key in hdu_noreprj[data_extension].header:
            if key in primary_hdr1:
                primary_hdr1[key] = hdu_noreprj[data_extension].header[key]

        # No change to second file header
        primary_hdr2 = hdu_noreprj[0].header
    else:
        # If the first file is not reprojected (resampled) then state it was only
        #   convolved in the primary header HISTORY
        primary_hdr1 = hdu_noreprj[0].header
        primary_hdr1['HISTORY'] = 'Convolution'
        primary_hdr1['HISTORY'] = '-- Map was convolved'
        primary_hdr1['HISTORY'] = '-- Parameter inputs'
        primary_hdr1['HISTORY'] = '  sigma_smooth = ' + str(sigma_smooth) + ' arcsec'
        primary_hdr1['HISTORY'] = '  sigma_beam = ' + str(sigma_beam) + ' arcsec'
        for key in kwargs:
            primary_hdr1['HISTORY'] = '  ' + key + ' = ' + str(kwargs[key])
        primary_hdr1['HISTORY'] = '--'
        primary_hdr1['HISTORY'] = ''

        # If the second file is reprojected (resampled) then state it was
        #   reprojected in the primary header HISTORY
        primary_hdr2 = hdu_reprj[0].header
        primary_hdr2['HISTORY'] = 'Reprojection'
        primary_hdr2['HISTORY'] = '-- Map was reprojected (resampled)'
        primary_hdr2['HISTORY'] = '-- Parameter inputs'
        primary_hdr2['HISTORY'] = '  reprojection_method = ' + str(reprojection_method)
        primary_hdr2['HISTORY'] = '  pixel_area_unit_type = ' + pixel_unit
        for key in kwargs:
            primary_hdr2['HISTORY'] = '  ' + key + ' = ' + str(kwargs[key])
        primary_hdr2['HISTORY'] = '--'
        primary_hdr2['HISTORY'] = ''
        # If the header had an extension name, remove it
        if 'EXTNAME' in primary_hdr2:
            primary_hdr2.remove('EXTNAME')
        # If a keyword is in both the primary and non-reprojected image header, replace
        #   the primary value with the non-reprojected value
        for key in hdu_noreprj[data_extension].header:
            if key in primary_hdr2:
                primary_hdr2[key] = hdu_noreprj[data_extension].header[key]

    # The first file which was convolved always gets the non-reprojected header
    convol_hd_data = hdu_noreprj[data_extension].header
    convol_hd_unc = hdu_noreprj[unc_extension].header

    # Generate HDUs for saving to output files
    if file_to_reproject == 1:
        primary1_hdu = fits.PrimaryHDU(header=primary_hdr1)
        data1_hdu = fits.ImageHDU(convol_data, header=convol_hd_data)
        unc1_hdu = fits.ImageHDU(convol_unc, header=convol_hd_unc)
        primary2_hdu = fits.PrimaryHDU(header=primary_hdr2)
        # If the second file was not reprojected (resampled) then it
        #   has no change from the input file
        data2_hdu = hdu_noreprj[data_extension]
        unc2_hdu = hdu_noreprj[unc_extension]
    else:
        primary1_hdu = fits.PrimaryHDU(header=primary_hdr1)
        data1_hdu = fits.ImageHDU(convol_data, header=convol_hd_data)
        unc1_hdu = fits.ImageHDU(convol_unc, header=convol_hd_unc)
        primary2_hdu = fits.PrimaryHDU(header=primary_hdr2)
        # If the second file was reprojected (resampled) then it
        #   gets the non-reprojected data and uncertainty's headers
        data2_hdu = fits.ImageHDU(data_reprj, header=hdu_noreprj[data_extension].header)
        unc2_hdu = fits.ImageHDU(unc_reprj, header=hdu_noreprj[unc_extension].header)

    # Save the reprojected and convolved images as extensions to the blank Primary
    hdul1_new = fits.HDUList([primary1_hdu, data1_hdu, unc1_hdu])
    hdul2_new = fits.HDUList([primary2_hdu, data2_hdu, unc2_hdu])
    hdul1_new.writeto(out_filename1, overwrite=True)
    hdul2_new.writeto(out_filename2, overwrite=True)
