def line_map(out_filename, filename, extensions, center_wavelength, velocity=0, revise_bounds=False, snr_limit=0,
             mcmc=False, **kwargs):
    """
    Wrapper function that reads a FITS file and fits an emission
      line with a Gaussian with the optional addition of up to a
      2nd degree polynomial. It then compiles the fits into a
      new FITS file containing the resulting line intensity,
      line intensity uncertainty, continuum, velocity, and
      FWHM maps.
    Parameters
    ----------
    out_filename : string
        A string containing the name of the resulting FITS file.
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
        A scalar containing the velocity of the object in km/s.
        If not specified a value of 0 is assumed.
    revise_bounds : boolean, optional
        A boolean that if set to True will refit the data using
        an initial fit's parameter ranges as new bounds for the
        parameter ranges.
    snr_limit : scalar, optional
        A scalar which is only used if 'revise_bounds' is True.
        It indicates a signal-to-noise level of the
        initial fit's line intensity below which
        data will not be considered when revising the bounds.
    mcmc : bool, optional
        A boolean specifying if an MCMC algorithm should be used to
        fit the model to the data. The MCMC algorithm uses the default
        emcee package (https://emcee.readthedocs.io/en/stable/user/install/).
        The initial state of the MCMC chain is the result from the non-linear
        least squares fit and the log-probability come from chisqr.
    kwargs
        Keyword arguments passed to the function line_fitting().
    """

    from spec_map_analysis.spectra_fitting import line_fitting
    from astropy.io import fits
    import numpy as np
    from spec_map_analysis.spectra_fitting import file_reader
    from copy import deepcopy

    # Read in the data and generate headers for output FITS files
    # Copy the kwargs dictionary and add in misc keywords for addition to the primary header HISTORY and ease
    #   of use in file_reader function
    kwargs_reader = deepcopy(kwargs)
    kwargs_reader['revise_bounds'] = revise_bounds
    kwargs_reader['snr_limit'] = snr_limit
    kwargs_reader['mcmc'] = mcmc
    fitting_data, primary_hdr, image_hdr = file_reader(filename, extensions, center_wavelength,
                                                       velocity=velocity, **kwargs_reader)

    # Fit the data, and if bounds are to be revised, do not fit with MCMC
    if revise_bounds:
        line_intensity, parameter = line_fitting(fitting_data, **kwargs)
    else:
        line_intensity, parameter = line_fitting(fitting_data, mcmc=mcmc, **kwargs)

    # If the keyword revise_bounds is set, refit the data using the current fit to further
    #   restrict the fitting bounds
    # Check if number of terms in the fit is specified. If not set to default of 3
    if 'nterms' in kwargs:
        nterms = kwargs['nterms']
    else:
        nterms = 3

    if revise_bounds:
        # Refit the data using the initial fits as better constraints on the Gaussian peak location and
        #   sigma ranges as to generate better fits.
        # Create bounds for each fit parameter based on previous high SNR fits. Exclude those with
        #   extremely high signal-to-noise as it is likely a artifact of the fitting
        snr = line_intensity['INTENSITY'] / line_intensity['UNCERTAINTY']
        snr[snr > 1e3] = 0
        snr_mask = snr > snr_limit

        vel = line_intensity['VELOCITY']
        width = line_intensity['FWHM']

        # Check if lower bounds were set. If set, use them for peak height, and continuum limits.
        #   Note: the revised bounds can only reduce the bound range from the input, and cannot expand it
        if 'lower_bounds' in kwargs:
            lower_bound = kwargs['lower_bounds']
            lower = np.array([lower_bound[0], np.nanmin(vel[snr_mask]), np.nanmin(width[snr_mask])])
            if nterms >= 4:
                lower = np.append(lower, lower_bound[3])
            if nterms >= 5:
                lower = np.append(lower, lower_bound[4])
            if nterms == 6:
                lower = np.append(lower, lower_bound[5])
        else:
            # If not set, use the input data and continuum to set peak and flat continuum. Do not restrict the
            #   linear or quadratic component of the continuum.
            lower = np.array([np.nanmin(fitting_data['FLUX']), np.nanmin(vel[snr_mask]), 0])
            if nterms >= 4:
                lower = np.append(lower, np.nanmin(line_intensity['CONTINUUM']))
            if nterms >= 5:
                lower = np.append(lower, -np.inf)
            if nterms == 6:
                lower = np.append(lower, -np.inf)

        # Check if upper bounds were set. If set, use them for peak height, and continuum limits.
        #   Note: the revised bounds can only reduce the bound range from the input, and cannot expand it
        if 'upper_bounds' in kwargs:
            upper_bound = kwargs['upper_bounds']
            upper = np.array([upper_bound[0], np.nanmax(vel[snr_mask]), np.nanmax(width[snr_mask])])
            if nterms >= 4:
                upper = np.append(upper, upper_bound[3])
            if nterms >= 5:
                upper = np.append(upper, upper_bound[4])
            if nterms == 6:
                upper = np.append(upper, upper_bound[5])
        else:
            # If not set, use the input data and continuum to set peak and flat continuum. Do not restrict the
            #   linear or quadratic component of the continuum.
            upper = np.array([np.nanmax(fitting_data['FLUX']), np.nanmax(vel[snr_mask]), np.nanmax(width[snr_mask])])
            if nterms >= 4:
                upper = np.append(upper, np.nanmax(line_intensity['CONTINUUM']))
            if nterms >= 5:
                upper = np.append(upper, np.inf)
            if nterms == 6:
                upper = np.append(upper, np.inf)

        # Refit the data using the new restricted bounds.
        kwargs['upper_bounds'] = upper
        kwargs['lower_bounds'] = lower
        line_intensity, parameter = line_fitting(fitting_data, mcmc=mcmc, **kwargs)

    # Construct the Primary extension
    primary_hdu = fits.PrimaryHDU(header=primary_hdr)

    # Construct the image extensions and assign appropriate units and extension names
    intensity = line_intensity['INTENSITY']
    image_hdr['BUNIT'] = 'W/m^2/pixel'
    image_hdr['EXTNAME'] = 'LINE_INTENSITY'
    intensity_hdu = fits.ImageHDU(intensity, header=image_hdr)

    uncertainty = line_intensity['UNCERTAINTY']
    image_hdr['BUNIT'] = 'W/m^2/pixel'
    image_hdr['EXTNAME'] = 'UNCERTAINTY'
    uncertainty_hdu = fits.ImageHDU(uncertainty, header=image_hdr)

    continuum = line_intensity['CONTINUUM']
    image_hdr['BUNIT'] = 'Jy/pixel'
    image_hdr['EXTNAME'] = 'CONTINUUM'
    continuum_hdu = fits.ImageHDU(continuum, header=image_hdr)

    velocity = line_intensity['VELOCITY']
    image_hdr['BUNIT'] = 'km/s/pixel'
    image_hdr['EXTNAME'] = 'VELOCITY'
    velocity_hdu = fits.ImageHDU(velocity, header=image_hdr)

    fwhm = line_intensity['FWHM']
    image_hdr['BUNIT'] = 'km/s/pixel'
    image_hdr['EXTNAME'] = 'FWHM'
    fwhm_hdu = fits.ImageHDU(fwhm, header=image_hdr)

    if mcmc:
        param = parameter['MCMC_CHAIN']
    else:
        param = parameter['PARAMETERS']
    param_hdr = intensity_hdu.header.copy()
    param_hdr['BUNIT'] = 'See comments'
    if mcmc:
        param_hdr['CUNIT4'] = 'parameter'
        param_hdr['CUNIT3'] = 'MCMC chain elements'
    else:
        param_hdr['CUNIT3'] = 'parameter'
    param_hdr['EXTNAME'] = 'PARAMETERS'
    param_hdr['COMMENT'] = 'The parameters in the axis containing the parameters are as follows:'
    param_hdr['COMMENT'] = 'First Index: Height of the Gaussian peak in Jy'
    param_hdr['COMMENT'] = 'Second Index: Frequency offset of the Gaussian peak from the rest-frame in Hz'
    param_hdr['COMMENT'] = 'Third Index: Sigma of the Gaussian in delta_frequency in Hz'
    if nterms >= 4:
        param_hdr['COMMENT'] = 'Fourth Index: 0th degree continuum component in Jy'
    if nterms >= 5:
        param_hdr['COMMENT'] = 'Fifth Index: 1st degree continuum component in Jy'
    if nterms == 6:
        param_hdr['COMMENT'] = 'Sixth Index: 2nd degree continuum component in Jy'
    param_hdu = fits.ImageHDU(param, header=param_hdr)

    if not mcmc:
        cov = parameter['COVARIANCE']
        cov_hdr = intensity_hdu.header.copy()
        cov_hdr['BUNIT'] = 'See comments'
        cov_hdr['CUNIT3'] = 'covariance'
        cov_hdr['CUNIT4'] = 'covariance'
        cov_hdr['EXTNAME'] = 'COVARIANCE'
        cov_hdr['COMMENT'] = 'The final two axes contain the covariance matrix of the fit parameters'
        cov_hdr['COMMENT'] = "See 'PARAMETERS' header for units of the parameters"
        cov_hdu = fits.ImageHDU(cov, header=cov_hdr)

    # Create a list of the primary, image, and table extensions and write to location
    if not mcmc:
        new_hdul = fits.HDUList([primary_hdu, intensity_hdu, uncertainty_hdu,
                                 continuum_hdu, velocity_hdu, fwhm_hdu, param_hdu, cov_hdu])
    else:
        new_hdul = fits.HDUList([primary_hdu, intensity_hdu, uncertainty_hdu,
                                 continuum_hdu, velocity_hdu, fwhm_hdu, param_hdu])
    new_hdul.writeto(out_filename, overwrite=True)
