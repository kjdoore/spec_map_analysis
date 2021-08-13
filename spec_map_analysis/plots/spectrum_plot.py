def spectrum_plot(out_filename, input_file, extensions, locations, location_unit='pixel',
                  center_wavelength=-1, velocity=0, fit_file=None, edge_buffer=None,
                  out_file_type='png', units=None, **kwargs):
    """
    Plots the spectra and uncertainty range for any amount of given
    pixels from a spectral data cube. Can also over-plot the fit to
    the spectra and include the residuals if a FITS file from the
    spectra_fitting pipeline is specified.
    Parameters
    ----------
    out_filename : string
        A string containing the base name of the file that the
        figure is to be saved. The string is appended to include
        the location and is saved as out_file_type.
    input_file : string
        A string containing the FITS file with the spectral data
        cube to be read.
    extensions : list of strings or integers
        A list of 3 or 4 string and/or integers containing the name
        or index of the extensions to be read. The order of the
        list must be 0) the flux data cube extension, 1) the flux
        error data cube extension, 2) the array of wavelengths
        extension, and optionally 3) the exposure map data cube
        extension. If the wavelength array is in a FITS table, a tuple
        can be given for the wavelength extension, which gives the
        table extension and table column name, respectively.
    locations : list of tuples
        A list of 2 element tuples containing the locations at which
        to plot the spectrum. Depending on the value of location_unit,
        the first element of the tuple must either be the Right Ascension
        of the location or the x-axis pixel location, and the second element
        must either be the Declination of the location or the y-axis pixel
        location.
    location_unit : string, optional
        A case-insensitive string containing the unit type for locations.
        The options are 'sky' and 'pixel' (default). 'sky' units refer
        to sky units as in an angle in reference to a given sky
        coordinate system, and 'pixel' units refer to the native pixel
        units of the map array. The unit will be the same for all locations.
    center_wavelength : scalar, optional
        A scalar containing the expected center wavelength of
        the spectral line as if observed in the rest-frame. Units
        of center_wavelength must match those specified in the units
        keyword or be in microns if not specified.
        If specified, a secondary x-axis will be plotted showing the
        velocity of the line.
    velocity : scalar, optional
        A scalar containing the velocity of the object in km/s
        to the chosen standard of rest.
        If not specified, a default value of 0 is used.
    fit_file : string, optional
         A string containing the FITS file with the fit to the data
         cube to be read. If not specified, then the fit to the
         spectrum will not be over-plotted.
    edge_buffer : two element list of integers, optional
        A two element list of integers that gives the number of
        wavelength indices to remove from the left (first element)
        and right (second element) edges of the wavelength range
        when performing fits. If not specified, no indices are removed
        from the wavelength range.
    out_file_type : string, optional
        A string containing the file type to save the image as. Default
        is 'png'. See https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.savefig.html
        for supported file types.
    units : list of strings, optional
        A two element list of strings containing the units to be included
        in the plots for the flux and wavelength. The first element is the
        flux unit, the second is the wavelength unit.
    kwargs
        Additional Keyword arguments passed to astropy.coordinates.SkyCoord
    """

    import matplotlib.pyplot as plt
    import numpy as np
    from spec_map_analysis.spectra_fitting import gauss_func as gf
    from astropy.io import fits
    from astropy import constants as const
    from spec_map_analysis.spectra_fitting import file_reader
    from astropy.coordinates import SkyCoord
    from astropy.wcs import WCS

    if units is None:
        units = ['Jy/pixel', '$\mu$m']

    # read in the spectral cube and place into convenient dictionary
    fitting_data, _, _ = file_reader(input_file, extensions, center_wavelength, velocity=velocity)

    # Compute the central wavelength of the line in the observed-frame
    wave_obs = np.sqrt((1 + fitting_data['VELOCITY'] / const.c.to('km/s').value) /
                       (1 - fitting_data['VELOCITY'] / const.c.to('km/s').value)) * fitting_data['CENTER_WAVE']

    # Place the data into individual arrays and convert the wavelength into velocity if wanted
    wave = fitting_data['WAVELENGTH']
    flux = fitting_data['FLUX']
    flux_unc = fitting_data['FLUX_UNC']
    if center_wavelength > 0:
        vel = ((wave / wave_obs) ** 2 - 1) / ((wave / wave_obs) ** 2 + 1) * const.c.to('km/s').value

    # Set flux to nan at the edges of the wavelength range to remove values from the plot
    if edge_buffer is not None:
        flux[0:edge_buffer[0], :, :] = float('nan')
        flux[-edge_buffer[1]::, :, :] = float('nan')
        flux_unc[0:edge_buffer[0], :, :] = float('nan')
        flux_unc[-edge_buffer[1]::, :, :] = float('nan')

    # If the fit data file is specified read it in and assign parameters to arrays. Also
    #   convert the wavelength into frequency since fits are performed in frequency space.
    if fit_file is not None:
        fit = fits.open(fit_file)
        param = fit['PARAMETERS'].data
        covar = fit['COVARIANCE'].data
        freq = const.c.to('um/s').value / fitting_data['WAVELENGTH']

    # Loop through and plot each spectrum.
    for i in range(0, len(locations)):
        # Generate output file name
        if location_unit.casefold() == 'pixel':
            output_file = out_filename + '_pixel_x' + str(locations[i][0]) + '_y' + str(locations[i][1])
        elif location_unit.casefold() == 'sky':
            output_file = out_filename + '_ra_' + str(locations[i][0]) + '_dec_' + str(locations[i][1])
        else:
            raise ValueError("location_unit must be 'pixel' or 'sky'")

        # Determine the pixel of the spectrum to plot if location given in sky units
        header = fits.getheader(input_file, extensions[0])
        if location_unit.casefold() == 'pixel':
            xpixel = locations[i][0]
            ypixel = locations[i][1]
        elif location_unit.casefold() == 'sky':
            sky = SkyCoord([locations[i]], **kwargs)
            # Since input_file is a cube, remove the spectral dimension
            header['NAXIS'] = 2
            for j in range(len(header.cards)-1, 0, -1):
                if '3' in header.cards[j][0]:
                    header.remove(header.cards[j][0])
            w = WCS(header)
            xpixel, ypixel = w.world_to_pixel(sky)
            xpixel = int(np.round(xpixel))
            ypixel = int(np.round(ypixel))
        else:
            raise ValueError("location_unit must be 'pixel' or 'sky'")

        # Only plot the data that has both data and uncertainty at a given wavelength. If not,
        #   then do not plot that wavelength.
        non_nan_flux = ~np.isnan(flux[:, ypixel, xpixel])
        non_nan_unc = ~np.isnan(flux_unc[:, ypixel, xpixel])
        non_nan = np.logical_and(non_nan_flux, non_nan_unc)

        # If the fit data file is given, select the parameters for the current pixel and generate its model
        #   spectrum, uncertainty, and the residual with the actual data.
        if fit_file is not None:
            par = param[:, ypixel, xpixel]
            cov = covar[:, :, ypixel, xpixel]

            mod_flux = gf.gauss_func(freq[non_nan], *par)
            flux_jacob = gf.gaussfunc_jacobian(freq[non_nan], *par)
            mod_flux_unc = np.diag(np.matmul(np.matmul(flux_jacob, cov), np.transpose(flux_jacob)))

            residual = (flux[non_nan, ypixel, xpixel] - mod_flux)

        # Generate the axes for the plot depending if the fit file was specified and
        #   residuals need to be plotted
        if fit_file is not None:
            fig, (ax0, ax2) = plt.subplots(2, 1, sharex='col')
            ax0.set_position([0.1, 0.3, 0.8, 0.6])
        else:
            fig, ax0 = plt.subplots(1, 1)

        # Plot the observed spectrum and give error range
        ax0.tick_params(direction='in')
        ax0.step(wave[non_nan], flux[non_nan, ypixel, xpixel], where='mid')
        ax0.fill_between(wave[non_nan], flux[non_nan, ypixel, xpixel] - flux_unc[non_nan, ypixel, xpixel],
                         flux[non_nan, ypixel, xpixel] + flux_unc[non_nan, ypixel, xpixel], alpha=0.4)
        ax0.plot([np.min(wave[non_nan]), np.max(wave[non_nan])], [0, 0], 'r--')
        ax0.set_ylabel('Flux ['+units[0]+']')

        # Add an upper x axis to be labeled with velocity
        if center_wavelength > 0:
            ax1 = ax0.twiny()

        # If the fit data file is given, plot the model and uncertainty
        if fit_file is not None:
            ax0.plot(wave[non_nan], mod_flux)
            ax0.fill_between(wave[non_nan], mod_flux - mod_flux_unc, mod_flux + mod_flux_unc, alpha=0.4)
            if center_wavelength > 0:
                ax1.set_position([0.1, 0.3, 0.8, 0.6])

        # Generate the velocity axis, but plot nothing. Just need the labels
        if center_wavelength > 0:
            ax1.tick_params(direction='in')
            ax1.set_xlabel('Velocity [km/s]')
            ax1.plot(vel[non_nan], flux[non_nan, ypixel, xpixel], ls=' ')

        # If the fit data file is given, plot the residuals
        if fit_file is not None:
            ax2.set_position([0.1, 0.1, 0.8, 0.2])
            ax2.tick_params(direction='in')
            ax2.step(wave[non_nan], residual, where='mid')
            ax2.fill_between(wave[non_nan], residual - flux_unc[non_nan, ypixel, xpixel],
                             residual + flux_unc[non_nan, ypixel, xpixel], alpha=0.4)
            ax2.plot([np.min(wave[non_nan]), np.max(wave[non_nan])], [0, 0], 'r--')
            ax2.set_xlabel('Wavelength ['+units[1]+']')
            ax2.set_ylabel('Residual')
        else:
            ax0.set_xlabel('Wavelength ['+units[1]+']')

        fig.savefig(output_file+'.'+out_file_type, dpi=500)
        plt.close('all')
