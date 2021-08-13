def line_fitting(fitting_data, nterms=3, exp_time_frac=0.0, edge_buffer=None,
                 upper_bounds=float('NaN'), lower_bounds=float('NaN'), mcmc=False, nwalkers=12,
                 nsteps=150, discard=50):
    """
    Fits spectral data cubes with a Gaussian with the optional
     addition of a polynomial up to degree 2. If nans are present,
     then they will be ignored in the fitting process. If not enough
     data is present to fit a spectrum, then no fit will be preformed.
    Parameters
    ----------
    fitting_data : dictionary
        A dictionary that contains the keys: 'WAVELENGTH', 'FLUX',
        'FLUX_UNC', 'VELOCITY', 'CENTER_WAVE', and optional 'EXPOSURE_MAP'.
        'WAVELENGTH' - array of length W containing the wavelengths
                       for each flux in microns
        'FLUX' - a (W, x, y) array containing the flux values at each
                 wavelength W and each pixel x and y in Jy
        'FLUX_UNC' - a (W, x, y) array containing the uncertainty on
                     the flux values at each wavelength W and each
                     pixel x and y in Jy
        'VELOCITY' - a scalar containing the velocity of the object
                     from the local standard of rest in units of km/s
        'CENTER_WAVE' - a scalar containing the center wavelength
                        of the line in microns that is to be fit as
                        if observed in the rest frame
        'EXPOSURE_MAP' (optional) - a (W, x, y) array containing the
                         exposure time at each each wavelength W and
                         each pixel x and y.
    nterms : integer, optional
        An integer value from 3 to 6 specifying the number of terms to
        include in the Gaussian function given by:
        f = a0*exp(-((x-a1)/a2)^2/2) + a3 + a4*x + a5*x^2
        If only a value of 3 is specified then only estimates for a
        Gaussian are returned (a0, a1, a2). If a value > 3 is
        specified then a polynomial of degree = nterms - 4 is added
        to the Gaussian up to a degree of 2 (a3, a4, a5). If a value
        is not specified, then a default value of 3 is used.
    exp_time_frac : scalar, optional
        A value between 0 and 1 that gives the fraction of the maximum
        exposure time below which data are not considered when
        performing fits. If not specified, the a default value of 0 is
        used.
    edge_buffer : two element list of integers, optional
        A two element list of integers that gives the number of
        wavelength indices to remove from the left (first element)
        and right (second element) edges of the wavelength range
        when performing fits. If not specified, no indices are removed
        from the wavelength range.
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
    mcmc : boolean, optional
        A boolean specifying if an MCMC algorithm should be used to
        fit the model to the data. The MCMC algorithm uses the default
        emcee package (https://emcee.readthedocs.io/en/stable/user/install/).
        The initial state of the MCMC chain is the result from the non-linear
        least squares fit and the log-probability come from chisqr.
    nwalkers : integer, optional
        The number of walkers in the MCMC ensemble. If mcmc is False
        this value is ignored
    nsteps : integer, optional
        The number of steps to run in the MCMC algorithm. If mcmc is False
        this value is ignored
    discard: integer, optional
        The beginning number of elements in the MCMC chain to discard
        as burn-in. If mcmc is False this value is ignored
    Returns
    -------
    intensity : dictionary
        A dictionary that contains the keys: 'INTENSITY',
        'UNCERTAINTY', 'CONTINUUM', 'VELOCITY', and 'FWHM'
        'INTENSITY' - a (x, y) array containing the line intensity
                      in W/m^2
        'UNCERTAINTY' - a (x, y) array containing the uncertainty on
                        the line intensity in W/m^2
        'CONTINUUM' - a (x, y) array containing the continuum flux
                      in Jy
        'VELOCITY' - a (x, y) array containing the velocity shift of
                     the line center from the galaxy rest-frame in km/s
        'FWHM' - a (x, y) array containing the the full width at half
                 maximum (FWHM) of the line in km/s
    parameter : dictionary
        A dictionary that depends upon the keyword mcmc:
        If mcmc=False then it contains the keys: 'FIT_PARAMETERS' and
        'COVARIANCE'.
        'PARAMETERS' - a (nterms, x, y) array containing the best fit
                       parameters for each pixel
        'COVARIANCE' - a (nterms, nterms, x, y) array containing the
                       covariance matrix from the best fit model for
                       each pixel
         If mcmc=True then it contains the key: 'MCMC_CHAIN'.
        'MCMC_CHAIN' - a (nterms, (nsteps-discard) * nwalkers, x, y) array
                       containing the MCMC chain for each parameters and pixel
   """

    import numpy as np
    from astropy import constants as const
    from spec_map_analysis.spectra_fitting import gauss_func as gf
    from scipy.optimize import curve_fit
    from spec_map_analysis.spectra_fitting import log_prob
    if mcmc:
        import emcee

    # Check to make sure inputs are of correct form
    # nterms needs to be an integer between 3 and 6. So, make integer and check range.
    nterms = int(nterms)
    if not 3 <= nterms <= 6:
        raise ValueError(
            ("nterms must be between 3 and 6; value given: {0}"
             ).format(nterms)
        )
    # exp_time_frac must be between 0 and 1
    if not 0 <= exp_time_frac <= 1:
        raise ValueError("exp_time_frac must be between 0 and 1")
    # edge_buffer values needs to be non-negative integers and combined are less than the
    #   number of wavelength values. So, convert to integer and check range.
    if edge_buffer is not None:
        edge_buffer[0] = int(edge_buffer[0])
        edge_buffer[1] = int(edge_buffer[1])
        if not 0 <= np.sum(edge_buffer) < fitting_data['WAVELENGTH'].size:
            raise ValueError(
                ("the sum of edge_buffer ({0}) must be less than the number of wavelength elements ({1})"
                 ).format(np.sum(edge_buffer), fitting_data['WAVELENGTH'].size)
            )
    # lower_bounds and upper_bounds must have nterms number of elements
    if np.sum(~np.isnan(lower_bounds)) != 0 and np.sum(~np.isnan(upper_bounds)) != 0:
        if not (np.sum(~np.isnan(lower_bounds)) == nterms and np.sum(~np.isnan(upper_bounds)) == nterms):
            raise ValueError(
                ("lower_bounds and upper_bounds must have {0} elements"
                 ).format(nterms)
            )
    # nsteps must be larger than discard
    if nsteps <= discard:
        raise ValueError("nsteps must be larger than discard")

    # Determine expected peak wavelength of observed line
    wave_obs = np.sqrt((1 + fitting_data['VELOCITY'] / const.c.to('km/s').value) /
                       (1 - fitting_data['VELOCITY'] / const.c.to('km/s').value)) * fitting_data['CENTER_WAVE']

    # Assign fitting data to arrays
    freq = const.c.to('um/s').value / fitting_data['WAVELENGTH']
    flux = fitting_data['FLUX']
    flux_unc = fitting_data['FLUX_UNC']

    # Set flux in each pixel to nan if the exposure time is given and less than a given fraction the max exposure time
    if 'EXPOSURE_MAP' in fitting_data:
        exposure_time = fitting_data['EXPOSURE_MAP']
        exp_time_mask = exposure_time <= np.max(exposure_time) * exp_time_frac
        flux[exp_time_mask] = float('nan')
        flux_unc[exp_time_mask] = float('nan')

    # Set flux to nan at the edges of the wavelength range to restrict noise from the edge of the detectors
    if edge_buffer is not None:
        flux[0:edge_buffer[0], :, :] = float('nan')
        flux[-edge_buffer[1]::, :, :] = float('nan')
        flux_unc[0:edge_buffer[0], :, :] = float('nan')
        flux_unc[-edge_buffer[1]::, :, :] = float('nan')


    # Create boolean saying if bounds are specified, and if so, set to appropriate units.
    if np.sum(~np.isnan(lower_bounds)) == nterms and np.sum(~np.isnan(upper_bounds)) == nterms:
        bounds_bool = True
        # Convert the bound peak location and FWHM input values to frequency units, and the FWHM to sigma
        #   FWHM = 2*sqrt(2*ln(2))*sigma
        # Note: the peak values are switched due to input in wavelength but calculation in frequency.
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
        lower_bounds[1], upper_bounds[1] = upper_bounds[1], lower_bounds[1]
    else:
        bounds_bool = False

    # Loop through each pixel and fit the spectra, then place the fit parameters and covariances into
    #   arrays with dimensions of [parameters, xpixel, ypixel] and [cov_matrix_dim1, cov_matrix_dim2, xpixel, ypixel]
    param_arr = np.empty((nterms, flux.shape[1], flux.shape[2]))
    cov_arr = np.empty((nterms, nterms, flux.shape[1], flux.shape[2]))
    # If the MCMC algorithm is to be used, create array to store each pixel's fit.
    if mcmc:
        mcmc_chain = np.empty((nterms, (nsteps-discard) * nwalkers, flux.shape[1], flux.shape[2]))

    for i in range(flux.shape[1]):
        for j in range(flux.shape[2]):

            # Index nan values that will be excluded from the fits
            non_nan_flux = ~np.isnan(flux[:, i, j])
            non_nan_unc = ~np.isnan(flux_unc[:, i, j])
            non_nan = np.logical_and(non_nan_flux, non_nan_unc)

            # Check to see if a pixel has enough spectral data for quality fit.
            #   If the amount of spectral data is less than the 1/3 of the total possible
            #   data points do not fit and assign fit parameters and uncertainties as nans
            # Also if there is no positive data do not fit and assign fit parameters and uncertainties
            #   as nans
            if np.sum(non_nan) < freq.size / 3 or np.nanmax(flux[:, i, j]) < 0:
                param_arr[:, i, j] = float('NaN')
                cov_arr[:, :, i, j] = float('NaN')
                if mcmc:
                    mcmc_chain[:, :, i, j] = float('NaN')
            else:
                # If no bounds are given, generate them. These bounds limit the parameters to ranges limited
                #   by the data itself (i.e., max peak height is slightly larger than the maximum flux value,
                #   velocity values cannot be outside of wavelength range, etc.)
                if not bounds_bool:
                    # Create bounds for each fit parameter for initial fit
                    lower_bounds, upper_bounds = gf.gaussfunc_bounds(freq[non_nan], flux[non_nan, i, j], nterms=nterms)

                # Generate initial guess for the parameters for fitting. Guesses are restricted to be
                #   within the specified bounds
                initial_guess = gf.gaussfunc_best_guess(freq[non_nan], flux[non_nan, i, j], nterms=nterms,
                                                        min_bounds=lower_bounds, max_bounds=upper_bounds)

                # Fit the data and derive the initial best fit parameters and the covariance matrix
                par, cov = curve_fit(gf.gauss_func, freq[non_nan], flux[non_nan, i, j], sigma=flux_unc[non_nan, i, j],
                                     jac=gf.gaussfunc_jacobian, p0=initial_guess, bounds=(lower_bounds, upper_bounds),
                                     maxfev=5000, absolute_sigma=True)
                param_arr[:, i, j] = par
                cov_arr[:, :, i, j] = cov

                # Fit using MCMC for better uncertainty estimation if keyword 'mcmc' is specified
                if mcmc:
                    # Create a Gaussian ball of starting positions for each ensemble walker
                    #   based off the maximum likelihood fit from the non-linear least squares
                    #   method. Scale the random value to each parameters order of magnitude.
                    mcmc_guess = np.tile(par, (nwalkers, 1)) + 1e-3 * np.random.randn(nwalkers, nterms) * \
                                 np.tile(10**np.floor(np.log10(np.abs(par))), (nwalkers, 1))
                    # Initialize the sampler
                    sampler = emcee.EnsembleSampler(nwalkers, nterms, log_prob,
                                                    args=[freq[non_nan], flux[non_nan, i, j], flux_unc[non_nan, i, j],
                                                          lower_bounds, upper_bounds])
                    # Run the sampler from the initial guesses for nsteps
                    state = sampler.run_mcmc(mcmc_guess, nsteps)
                    # Compress each walker into one chain, excluding each walker's discard number
                    #   of first steps as a burn-in phase
                    mcmc_chain[:, :, i, j] = np.transpose(sampler.get_chain(flat=True, discard=discard))

    # Place the parameters and covariance or MCMC chain into a dictionary
    if mcmc:
        parameter = {'MCMC_CHAIN': mcmc_chain}
    else:
        parameter = {'PARAMETERS': param_arr, 'COVARIANCE': cov_arr}

    # Calculate the line intensity and uncertainty.
    # Line intensity of the Gaussian function (a0*exp(-1/2*((nu-a1)/a2)^2) is calculated as line_int=a0*a2*sqrt(2*pi)
    #   This is derived from integral(exp(-a*(x+b)^2),-inf,inf) = sqrt(pi/a), where a=1/(2*a2^2) from the fitted
    #   Gaussian function. Propagating uncertainty with this equation gives:
    #   sig_line_int=sqrt(2*pi*a2^2*sig_a0^2+2*pi*a0^2*sig_a2^2+2*pi*a0*a2*cov_a0_a2)
    # Resulting units are in W m^-2
    l_int = param_arr[0, :, :] * param_arr[2, :, :] * np.sqrt(2 * np.pi) * 1e-26
    l_int_unc = np.sqrt((2 * np.pi * param_arr[2, :, :] ** 2 * cov_arr[0, 0, :, :] +
                         2 * np.pi * param_arr[0, :, :] ** 2 * cov_arr[2, 2, :, :] +
                         4 * np.pi * param_arr[0, :, :] * param_arr[2, :, :] * cov_arr[0, 2, :, :]) * 1e-26 * 1e-26)
    if mcmc:
        l_int_mcmc = np.median(mcmc_chain[0, :, :, :] * mcmc_chain[2, :, :, :] * np.sqrt(2 * np.pi) * 1e-26, axis=0)
        l_int_unc_mcmc = np.std(mcmc_chain[0, :, :, :] * mcmc_chain[2, :, :, :] * np.sqrt(2 * np.pi) * 1e-26, axis=0)

    # Calculate the continuum at the peak location for any value of nterms. It will be 0 for nterms = 3,
    #   and will have the flux density value at the peak location for nterms > 3
    if nterms <= 3:
        continuum = np.zeros((flux.shape[1], flux.shape[2]))
    else:
        continuum = np.polynomial.polynomial.polyval(param_arr[1, :, :], param_arr[3::, :, :], tensor=False)
    if mcmc:
        if nterms <= 3:
            continuum_mcmc = np.zeros((flux.shape[1], flux.shape[2]))
        else:
            continuum_mcmc = np.median(np.polynomial.polynomial.polyval(mcmc_chain[1, :, :, :],
                                                                        mcmc_chain[3::, :, :, :], tensor=False), axis=0)

    # Calculate the offset of the line peak from the expected wavelength in terms of velocity.
    #   The peak location in the param_arr variable is in frequency. Convert to wavelength and
    #   determine the velocity via the relativistic Doppler shift.
    peak_wave = const.c.to('um/s').value / param_arr[1, :, :]
    velocity = ((peak_wave / wave_obs) ** 2 - 1) / ((peak_wave / wave_obs) ** 2 + 1) * const.c.to('km/s').value
    if mcmc:
        peak_wave = const.c.to('um/s').value / mcmc_chain[1, :, :, :]
        velocity_mcmc = np.median(((peak_wave / wave_obs) ** 2 - 1) / ((peak_wave / wave_obs) ** 2 + 1) *
                                  const.c.to('km/s').value, axis=0)

    # Calculate the full width at half maximum (FWHM) of the line. Fit parameter is the sigma of the Gaussian.
    #   So the FWHM is 2*sqrt(2*ln(2))*sigma. Since a2 = sigma, FWHM = 2*sqrt(2*ln(2))*a2. Units of a2
    #   are delta_nu, so convert to km/s via delta_nu/nu_peak = delta_v/c
    fwhm = 2 * np.sqrt(2 * np.log(2)) * param_arr[2, :, :]
    fwhm = fwhm / param_arr[1, :, :] * const.c.to('km/s').value
    if mcmc:
        fwhm_mcmc = 2 * np.sqrt(2 * np.log(2)) * mcmc_chain[2, :, :, :]
        fwhm_mcmc = np.median(fwhm_mcmc / mcmc_chain[1, :, :, :] * const.c.to('km/s').value, axis=0)

    # Place values into dictionary
    if mcmc:
        line_intensity = {'INTENSITY': l_int_mcmc, 'UNCERTAINTY': l_int_unc_mcmc, 'CONTINUUM': continuum_mcmc,
                          'VELOCITY': velocity_mcmc, 'FWHM': fwhm_mcmc}
        return line_intensity, parameter
    else:
        line_intensity = {'INTENSITY': l_int, 'UNCERTAINTY': l_int_unc, 'CONTINUUM': continuum,
                          'VELOCITY': velocity, 'FWHM': fwhm}
        return line_intensity, parameter
