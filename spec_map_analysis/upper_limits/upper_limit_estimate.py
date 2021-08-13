def upper_limit_estimate(filename, extensions, center_wavelength, velocity, max_height, peak_velocity,
                         fwhm, mask=None, nsim=100, edge_buffer=0, exp_time_frac=0,
                         nterms=3, lower_bounds=float('NaN'), upper_bounds=float('NaN'), **kwargs):
    """
    Reads a FITS file with specified extensions, and returns
      an estimate on the upper limit of the line intensity
      required for a 1-sigma detection. Estimates are made
      when simulated sources placed on the image result in
      detections comparable to the simulated source's true
      line intensity.
    Parameters
    ----------
    filename : string
        A string containing the FITS file that upper limits
        are to be estimated.
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
    velocity : scalar
        A scalar containing the velocity of the object in km/s
        to the chosen standard of rest.
    max_height : scalar
        A scalar containing the maximum peak height of the
        simulated Gaussian in units of Jy. Larger values
        are required for noisier data.
    peak_velocity : scalar
        A scalar containing the velocity offset of the
        simulated line in km/s
    fwhm : scalar
        A scalar containing the full width at half maximum
        (fwhm) of the simulated line in km/s
    mask : boolean array
        A two dimensional boolean array with the same shape as the input
        file data. A True value indicates the corresponding element
        of the input array are to not have a simulated source inserted in
        that element.
    nsim : integer, optional
        An integer containing the number of simulated sources to
        individually input into the input image. Not all simulated sources
        result in an upper limit estimate if the chosen location already
        has a detectable source.
    edge_buffer : integer, optional
        An integer value that gives the number of wavelength indices to
        remove from both edges of the wavelength range when performing
        fits. If not specified, the a default value of 0 is used.
    exp_time_frac : scalar, optional
        A scalar between 0 and 1 that gives the fraction of the maximum
        exposure time below which data are not considered when
        performing fits. If not specified, the a default value of 0 is
        used.
    nterms : integer, optional
        An integer value from 3 to 6 specifying the number of terms to
        include in the Gaussian function given by:
        f = a0*exp(-((x-a1)/a2)^2/2) + a3 + a4*x + a5*x^2
        If only a value of 3 is specified then only estimates for a
        Gaussian are returned (a0, a1, a2). If a value > 3 is
        specified then a polynomial of degree = nterms - 4 is added
        to the Gaussian up to a degree of 2 (a3, a4, a5). If a value
        is not specified, then a default value of 3 is used.
    lower_bounds : 1-D array-like, optional
        An array of length nterms giving the lower bounds for the
        non-linear least squares fitting algorithm, where the potential
        values are (1) the height of the Gaussian in Jy, (2) the velocity
        offset of the peak in km/s, (3) the FWHM of the Gaussian in km/s,
        (4) the continuum in Jy, (5) the slope of the continuum, and (6)
        the quadratic continuum coefficient (it is recommended to set
        this to 0 if nterms=6 is used)
    upper_bounds : 1-D array-like, optional
        An array of length nterms giving the upper bounds for the
        non-linear least squares fitting algorithm, where the potential
        values are (1) the height of the Gaussian in Jy, (2) the velocity
        offset of the peak in km/s, (3) the FWHM of the Gaussian in km/s,
        (4) the continuum in Jy, (5) the slope of the continuum, and (6)
        the quadratic continuum coefficient (it is recommended to set
        this to 0 if nterms=6 is used)
    kwargs
        Keyword arguments passed to gaussfunc_bounds
    Returns
    -------
    upper_limit : array
        An array containing a distribution of upper limit
        estimates for the image.
    """

    from spec_map_analysis.upper_limits import simulated_source
    import numpy as np
    from astropy import constants as const
    from spec_map_analysis.spectra_fitting import gauss_func as gf
    from scipy.optimize import curve_fit
    from spec_map_analysis.spectra_fitting import file_reader

    # Read in the data using the extension names
    fitting_data, _, _ = file_reader(filename, extensions, center_wavelength, velocity=velocity)

    # Determine the observed-frame wavelength of the line
    wave_obs = np.sqrt((1 + fitting_data['VELOCITY'] / const.c.to('km/s').value) /
                       (1 - fitting_data['VELOCITY'] / const.c.to('km/s').value)) * fitting_data['CENTER_WAVE']

    # Assign fitting data to arrays
    freq = const.c.to('um/s').value / fitting_data['WAVELENGTH']
    flux = fitting_data['FLUX']
    flux_unc = fitting_data['FLUX_UNC']

    # If exposure map is specified, then limit the fitted area to only be only spectral and spatial
    #   pixels with more exposure time than the given fraction of the max exposure time
    if 'EXPOSURE_MAP' in fitting_data:
        exposure_time = fitting_data['EXPOSURE_MAP']
        exp_time_mask = exposure_time <= np.max(exposure_time) * exp_time_frac
        flux[exp_time_mask] = float('nan')
        flux_unc[exp_time_mask] = float('nan')

    # Set flux to nan at the edges of the wavelength range to restrict noise from the edge of the detectors
    if edge_buffer != 0:
        flux[0:edge_buffer[0], :, :] = float('nan')
        flux[-edge_buffer[1]::, :, :] = float('nan')
        flux_unc[0:edge_buffer[0], :, :] = float('nan')
        flux_unc[-edge_buffer[1]::, :, :] = float('nan')

    # Set flux to nan where the mask is True
    if mask is not None:
        flux[:, mask] = float('nan')
        flux_unc[:, mask] = float('nan')

    # Create boolean saying if bounds are specified, and if so, set to appropriate units.
    if np.sum(~np.isnan(lower_bounds)) == nterms and np.sum(~np.isnan(upper_bounds)) == nterms:
        bounds_bool = True
        # Convert the bound peak location and FWHM input values to frequency units, and the FWHM to sigma
        #   FWHM = 2*sqrt(2*ln(2))*sigma
        lower_bounds[1] = const.c.to('um/s').value / (np.sqrt((1 + lower_bounds[1] / const.c.to('km/s').value) /
                                                              (1 - lower_bounds[1] / const.c.to('km/s').value))
                                                      * wave_obs)
        lower_bounds[2] = lower_bounds[2] / (2 * np.sqrt(2 * np.log(2))) / \
                          const.c.to('km/s').value * lower_bounds[1]
        upper_bounds[1] = const.c.to('um/s').value / (np.sqrt((1 + upper_bounds[1] / const.c.to('km/s').value) /
                                                              (1 - upper_bounds[1] / const.c.to('km/s').value))
                                                      * wave_obs)
        upper_bounds[2] = upper_bounds[2] / (2 * np.sqrt(2 * np.log(2))) / \
                          const.c.to('km/s').value * upper_bounds[1]
        # Note: the peak values are switched due to input in wavelength but calculation in frequency.
        lower_bounds[1], upper_bounds[1] = upper_bounds[1], lower_bounds[1]
    else:
        bounds_bool = False

    # Generate arrays for later updating
    height = np.arange(0, max_height, 0.01)
    upper_limit = np.zeros(nsim)

    # Generate array where simulated source at a pixel has enough
    #   spectral data for quality fit. Only select if the amount of spectral
    #   data is more than the 1/3 of the total possible data points and
    #   there is positive data.
    non_nan_flux = ~np.isnan(flux)
    non_nan_unc = ~np.isnan(flux_unc)
    non_nan = np.logical_and(non_nan_flux, non_nan_unc)
    simulatable_locations = np.logical_and(np.sum(non_nan * 1, axis=0) > freq.size / 3, np.nanmax(flux, axis=0) > 0)

    # Loop for number of simulated source to include
    for j in range(0, upper_limit.size):
        l_int = np.zeros(height.size)
        l_int_unc = np.zeros(height.size)
        l_int_diff = np.zeros(height.size)
        t_int = np.zeros(height.size)

        # Randomly select allowable simulated source location
        idx = np.flatnonzero(simulatable_locations)
        sim_loc = np.unravel_index(np.random.choice(idx, 1), (flux.shape[1], flux.shape[2]))

        # Loop through various simulated peak heights to find the upper limit at the given source location
        for i in range(1, height.size):
            # Insert simulated source
            flux_sim, true_intensity = simulated_source(fitting_data['WAVELENGTH'], flux, center_wavelength,
                                                        sim_loc, [1, 1], height[i], peak_velocity, fwhm,
                                                        velocity=velocity)

            # Copy the spectra at the simulated source location for fitting
            flux_temp = flux_sim[:, sim_loc[0], sim_loc[1]].flatten()
            flux_unc_temp = flux_unc[:, sim_loc[0], sim_loc[1]].flatten()
            non_nan_flux_temp = ~np.isnan(flux_temp)
            non_nan_unc_temp = ~np.isnan(flux_unc_temp)
            non_nan_temp = np.logical_and(non_nan_flux_temp, non_nan_unc_temp)

            # If no bounds are given, generate them.
            if not bounds_bool:
                # Create bounds for each fit parameter for initial fit
                lower_bounds, upper_bounds = gf.gaussfunc_bounds(freq[non_nan_temp], flux_temp[non_nan_temp],
                                                                 nterms=nterms)
            # Generate initial guess for the parameters for fitting.
            initial_guess = gf.gaussfunc_best_guess(freq[non_nan_temp], flux_temp[non_nan_temp], nterms=nterms,
                                                    min_bounds=lower_bounds, max_bounds=upper_bounds)

            # Fit the data and derive the best fit parameters and the covariance matrix
            par, cov = curve_fit(gf.gauss_func, freq[non_nan_temp], flux_temp[non_nan_temp],
                                 sigma=flux_unc_temp[non_nan_temp], jac=gf.gaussfunc_jacobian, p0=initial_guess,
                                 bounds=(lower_bounds, upper_bounds), maxfev=5000, absolute_sigma=True)

            # Calculate the true line intensity, and measured line intensity and uncertainty
            t_int[i] = true_intensity[sim_loc[0], sim_loc[1]]
            l_int[i] = par[0] * par[2] * np.sqrt(2 * np.pi) * 1e-26
            l_int_unc[i] = np.sqrt((2 * np.pi * par[2] ** 2 * cov[0, 0] +
                                    2 * np.pi * par[0] ** 2 * cov[2, 2] +
                                    4 * np.pi * par[0] * par[2] * cov[0, 2]) * 1e-26 * 1e-26)
            # Calculate the ratio of the previous simulation to the new simulation. This will
            #   be used to determine when the line is detected.
            l_int_diff[i] = l_int[i]/l_int[i-1]

        # Check for where simulated source is detected
        # A source is considered detected when the SNR > 3 (Low SNR < 1 means
        #   the simulated source is not detected from the background) and the
        #   ratio between the previous simulated source's height spikes (This
        #   spike iss formed when the algorithm finally detects the source).
        snr = l_int/l_int_unc
        upper_limit_loc = np.where(np.logical_and(snr > 3, l_int_diff > 3))
        if upper_limit_loc[0].size != 0:
            upper_limit[j] = l_int[np.min(upper_limit_loc)]

    # Remove upper limit estimates that are 0, which occurs when
    #   a source is detected when no simulated source has been added
    #   (i.e., there is already a detectable line at the location)
    upper_limit = upper_limit[upper_limit != 0]

    return upper_limit
