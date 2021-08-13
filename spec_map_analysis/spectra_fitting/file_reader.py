def file_reader(filename, extensions, center_wavelength, velocity=0, **kwargs):
    """
    Reads a FITS file with specified extensions, compiles
      the extensions to be used for fitting the spectra, and
      generates new headers for the line intensity maps.
    Parameters
    ----------
    filename : string
        A string containing the FITS file to be read.
    extensions : list of strings or integers
        A list of 3 or 4 string and/or integers containing the name
        or index of the extensions to be read. The order of the
        list must be 0) the flux data cube extension, 1) the flux
        error data cube extension, 2) the array of wavelengths
        extension, and optionally 3) the exposure map data cube
        extension. If the wavelength array is in a FITS table, a tuple
        can be given for the wavelength extension, which gives the
        table extension and table column name, respectively.
    center_wavelength : scalar
        A scalar containing the center wavelength of the line
        in microns that is to be fit as if observed in the
        rest frame
    velocity : scalar, optional
        A scalar containing the velocity of the object in km/s
        to the chosen standard of rest.
        If not specified, a default value of 0 is used.
    kwargs
        Keyword arguments passed to the fitting function that are
        added to the HISTORY in the Primary header.
    Returns
    -------
    fitting_data : dictionary
        A dictionary that contains the keys: 'WAVELENGTH', 'FLUX',
        'FLUX_UNC', 'VELOCITY', 'CENTER_WAVE', and optional 'EXPOSURE_MAP'.
        'WAVELENGTH' - array of length W containing the wavelengths
                       for each flux
        'FLUX' - a (W, x, y) array containing the flux values at each
                 wavelength W and each pixel x and y.
        'FLUX_UNC' - a (W, x, y) array containing the uncertainty on
                     the flux values at each wavelength W and each
                     pixel x and y.
        'VELOCITY' - a single valued array containing the velocity of
                     the object from the local standard of rest in
                     units of km/s
        'CENTER_WAVE' - a single valued array containing the center
                        wavelength of the line as if observed in the
                        rest frame to be fit
        'EXPOSURE_MAP' (optional) - a (W, x, y) array containing the
                         exposure time at each each wavelength W and
                         each pixel x and y.
    primary_hdr : Primary FITS header
        The Primary FITS header to be attached to the final FITS file
        after the fitting of the spectra
    image_hdr : Image extension FITS header
        The Image extensions' FITS header to be attached to each image
        extension of final FITS file
    """

    from astropy.io import fits

    # Read the file
    hdul = fits.open(filename)
    # Read the Primary header and an image header
    primary_hdr = hdul[0].header
    image_header = hdul[1].header

    # Convert file data to usable arrays
    nextens = len(extensions)
    flux = hdul[extensions[0]].data
    flux_unc = hdul[extensions[1]].data
    if type(extensions[2]) == tuple:
        wavelength = hdul[extensions[2][0]].data[extensions[2][1]]
    else:
        wavelength = hdul[extensions[2]].data
    if nextens == 4:
        exposure_time = hdul[extensions[3]].data

    # Put all fitting data into convenient dictionary for ease of use
    fitting_data = {'WAVELENGTH': wavelength, 'FLUX': flux, 'FLUX_UNC': flux_unc,
                    'VELOCITY': velocity, 'CENTER_WAVE': center_wavelength}
    if nextens == 4:
        fitting_data['EXPOSURE_MAP'] = exposure_time

    # Create the headers for final data products. The Primary header is the same as the
    #   input Primary header except it does not have the third spectral dimension, and it
    #   has the added HISTORY comments.
    if 'EXTNAME' in primary_hdr:
        primary_hdr.remove('EXTNAME')
    if 'CTYPE3' in primary_hdr:
        primary_hdr.remove('CTYPE3')
    if 'CUNIT3' in primary_hdr:
        primary_hdr.remove('CUNIT3')
    if 'CRPIX3' in primary_hdr:
        primary_hdr.remove('CRPIX3')
    if 'CRVAL3' in primary_hdr:
        primary_hdr.remove('CRVAL3')
    if 'CDELT3' in primary_hdr:
        primary_hdr.remove('CDELT3')
    primary_hdr['HISTORY'] = 'Line Intensity Calculation'
    primary_hdr['HISTORY'] = '-- Lines were fit with a Gaussian + continuum'
    primary_hdr['HISTORY'] = '-- Parameter inputs for line fitting'
    primary_hdr['HISTORY'] = '  velocity = ' + str(velocity) + ' km/s'
    primary_hdr['HISTORY'] = '  center_wavelength = ' + str(center_wavelength) + ' um'
    # add all additional keywords used for fitting to the HISTORY
    for key in kwargs:
        primary_hdr['HISTORY'] = '  ' + key + ' = ' + str(kwargs[key])
    primary_hdr['HISTORY'] = '--'
    primary_hdr['HISTORY'] = ''

    # Construct a Image header that will be uniform for each image extension
    #   and across different instruments
    # Generate blank standard image extension header
    image_hdr = fits.ImageHDU().header

    # Set the units to blank, will be updated for each extension later
    image_hdr['BUNIT'] = ('', 'Data units')

    # If the Equinox units are in the primary header, set to match. If not check an
    #   image header and set to match. If neither, assume J2000 and print warning.
    if 'EQUINOX' in primary_hdr:
        image_hdr['EQUINOX'] = (primary_hdr['EQUINOX'], 'Equinox of celestial coordinate system')
    elif 'EQUINOX' in image_header:
        image_hdr['EQUINOX'] = (image_header['EQUINOX'], 'Equinox of celestial coordinate system')
    else:
        image_hdr['EQUINOX'] = (2000.0, 'Equinox of celestial coordinate system')
        print('WARNING: No equinox of celestial coordinate system found in file headers')

    # Read required axis header info from primary or image headers and insert into
    #   final Image header
    if 'CTYPE1' in primary_hdr:
        image_hdr['CTYPE1'] = (primary_hdr['CTYPE1'], 'Axis 1 type and projection')
    elif 'CTYPE1' in image_header:
        image_hdr['CTYPE1'] = (image_header['CTYPE1'], 'Axis 1 type and projection')
    else:
        print('WARNING: CTYPE1 not found in file headers')
    if 'CTYPE2' in primary_hdr:
        image_hdr['CTYPE2'] = (primary_hdr['CTYPE2'], 'Axis 2 type and projection')
    elif 'CTYPE2' in image_header:
        image_hdr['CTYPE2'] = (image_header['CTYPE2'], 'Axis 2 type and projection')
    else:
        print('WARNING: CTYPE2 not found in file headers')

    if 'CUNIT1' in primary_hdr:
        image_hdr['CUNIT1'] = (primary_hdr['CUNIT1'], 'Axis 1 units')
    elif 'CUNIT1' in image_header:
        image_hdr['CUNIT1'] = (image_header['CUNIT1'], 'Axis 1 units')
    else:
        print('WARNING: CUNIT1 not found in file headers')
    if 'CUNIT2' in primary_hdr:
        image_hdr['CUNIT2'] = (primary_hdr['CUNIT2'], 'Axis 2 units')
    elif 'CUNIT2' in image_header:
        image_hdr['CUNIT2'] = (image_header['CUNIT2'], 'Axis 2 units')
    else:
        print('WARNING: CUNIT2 not found in file headers')

    if 'CRPIX1' in primary_hdr:
        image_hdr['CRPIX1'] = (primary_hdr['CRPIX1'], 'Axis 1 reference pixel position')
    elif 'CRPIX1' in image_header:
        image_hdr['CRPIX1'] = (image_header['CRPIX1'], 'Axis 1 reference pixel position')
    else:
        print('WARNING: CRPIX1 not found in file headers')
    if 'CRPIX2' in primary_hdr:
        image_hdr['CRPIX2'] = (primary_hdr['CRPIX2'], 'Axis 2 reference pixel position')
    elif 'CRPIX2' in image_header:
        image_hdr['CRPIX2'] = (image_header['CRPIX2'], 'Axis 2 reference pixel position')
    else:
        print('WARNING: CRPIX2 not found in file headers')

    if 'CRVAL1' in primary_hdr:
        image_hdr['CRVAL1'] = (primary_hdr['CRVAL1'], 'CTYPE1 in CUNIT1 at reference pixel')
    elif 'CRVAL1' in image_header:
        image_hdr['CRVAL1'] = (image_header['CRVAL1'], 'CTYPE1 in CUNIT1 at reference pixel')
    else:
        print('WARNING: CRVAL1 not found in file headers')
    if 'CRVAL2' in primary_hdr:
        image_hdr['CRVAL2'] = (primary_hdr['CRVAL2'], 'CTYPE2 in CUNIT2 at reference pixel')
    elif 'CRVAL2' in image_header:
        image_hdr['CRVAL2'] = (image_header['CRVAL2'], 'CTYPE2 in CUNIT2 at reference pixel')
    else:
        print('WARNING: CRVAL2 not found in file headers')

    if 'CDELT1' in primary_hdr:
        image_hdr['CDELT1'] = (primary_hdr['CDELT1'], 'Axis 1 pixel scale')
    elif 'CDELT1' in image_header:
        image_hdr['CDELT1'] = (image_header['CDELT1'], 'Axis 1 pixel scale')
    else:
        print('WARNING: CDELT1 not found in file headers')
    if 'CDELT2' in primary_hdr:
        image_hdr['CDELT2'] = (primary_hdr['CDELT2'], 'Axis 2 pixel scale')
    elif 'CDELT2' in image_header:
        image_hdr['CDELT2'] = (image_header['CDELT2'], 'Axis 2 pixel scale')
    else:
        print('WARNING: CDELT2 not found in file headers')

    if 'CROTA2' in primary_hdr:
        image_hdr['CROTA2'] = (primary_hdr['CROTA2'], 'Rotation angle')
    elif 'CROTA2' in image_header:
        image_hdr['CROTA2'] = (image_header['CROTA2'], 'Rotation angle')
    else:
        image_hdr['CROTA2'] = (-0.0, 'Rotation angle')
        print('WARNING: CROTA2 not found in file headers, assuming 0')

    # Set the extension name to blank, will be updated for each extension later
    image_hdr['EXTNAME'] = ('', 'Extension name')

    return fitting_data, primary_hdr, image_hdr
