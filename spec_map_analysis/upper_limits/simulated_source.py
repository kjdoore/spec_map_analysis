def simulated_source(wavelength, flux, center_wavelength, sim_loc, sigma,
                     height, peak_velocity, fwhm, velocity=0):
    """
    Simulates a Gaussian line source both spatially and spectrally
      for the given parameters and adds it to the input flux data
      cube.
    Parameters
    ----------
    wavelength : 1-D array_like
        An array of length W containing the wavelengths
        for each flux in microns
    flux : 3-D array_like
        A (W, x, y) array containing the flux values at each
        wavelength W and each pixel x and y in Jy
    center_wavelength : scalar
        A scalar containing the center wavelength of the simulated
        line in microns that is to be added as if observed in the
        rest frame
    sim_loc : array_like
        A two element array containing the spatial location
        in pixels where the simulated source is to be inserted.
        The first element is the x location, and the second is
        y location.
    sigma : array_like
        A two element array containing the spatial standard
        deviation (sigma) in pixels of the simulated Gaussian
        source. The first element is the x sigma, and the
        second is y sigma. There is no option for covariance,
        and the source is always circular.
    height : scalar
        A scalar containing the peak height of the
        simulated Gaussian in units of Jy.
    peak_velocity : scalar
        A scalar containing the velocity offset of the
        simulated line in km/s
    fwhm : scalar
        A scalar containing the full width at half maximum
        (fwhm) of the simulated line in km/s
    velocity : scalar, optional
        A scalar containing the velocity of the object in km/s.
        If not specified a value of 0 is assumed.
    Returns
    -------
    flux : array
        a (W, x, y) array containing the flux with the added
        simulated source values at each wavelength W and each
        pixel x and y in Jy
    true_intensity : array
        a (x, y) array containing the true line intensity of the
        simulated source at each pixel x and y in W/s^2
    """

    import numpy as np
    from astropy import constants as const
    from spec_map_analysis.spectra_fitting import gauss_func as gf

    freq = const.c.to('um/s').value / wavelength

    # Generate grids of x and y values at each pixel
    x_array = np.transpose(np.tile(np.arange(0, flux.shape[1]), (flux.shape[2], 1))) - sim_loc[0]
    y_array = np.tile(np.arange(0, flux.shape[2]), (flux.shape[1], 1)) - sim_loc[1]

    # Compute the spatial flux normalized amplitude
    spatial_flux_amp = np.exp(-0.5 * ((x_array/sigma[0])**2 + (y_array/sigma[1])**2))

    # Determine the observed frame line peak wavelength and use it to determine frequency
    #   at the simulated sources peak
    wave_obs = np.sqrt((1 + velocity / const.c.to('km/s').value) /
                       (1 - velocity / const.c.to('km/s').value)) * center_wavelength
    peak = const.c.to('um/s').value / (np.sqrt((1 + peak_velocity/const.c.to('km/s').value)/
                                               (1 - peak_velocity/const.c.to('km/s').value)) * wave_obs)

    # Convert the FWHM in velocity to sigma in delta_frequency
    fwhm = fwhm / (2 * np.sqrt(2 * np.log(2))) / const.c.to('km/s').value * peak

    # Add the new source to the existing data cube.
    flux_sim = flux + 0
    for i in range(0, flux.shape[1]):
        for j in range(0, flux.shape[2]):
            flux_sim[:, i, j] = flux[:, i, j] + gf.gauss_func(freq, height*spatial_flux_amp[i, j], peak, fwhm)

    # Create a map containing the true intensity of the simulated source at each pixel in W/m^2
    true_intensity = height*spatial_flux_amp * fwhm * np.sqrt(2 * np.pi) * 1e-26

    return flux_sim, true_intensity
