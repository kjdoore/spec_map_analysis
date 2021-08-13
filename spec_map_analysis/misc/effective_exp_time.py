def effective_exp_time(out_filename, filename, extensions, fit_file):
    """
    A function that reads in a FITS file containing a
      spectral map's exposure map and a FITS file containing
      the fit parameters from the fit spectral map. The
      effective exposure time for each pixel is then computed
      based on the mean exposure time of the spectral pixels
      that contain the FWHM of the fit Gaussian.
    Parameters
    ----------
    out_filename : string
        A string containing the output file that contains
        the effective exposure time map.
    filename : string
        A string containing the base FITS file to be read.
    extensions : list of strings or integers
        A list of 4 string and/or integers containing the name
        or index of the extensions to be read. The order of the
        list must be 0) the flux data cube extension, 1) the flux
        error data cube extension, 2) the array of wavelengths
        extension, and 3) the exposure map data cube extension.
        If the wavelength array is in a FITS table, a tuple
        can be given for the wavelength extension, which gives the
        table extension and table column name, respectively.
    fit_file : string, optional
        A string containing the FITS file with the fit to the data
        cube to be read.
    """

    from astropy.io import fits
    from spec_map_analysis.spectra_fitting import file_reader
    from astropy import constants as const
    import numpy as np

    # Read in the original data
    fitting_data, primary_hdr, image_hdr = file_reader(filename, extensions, None)

    # Read in the parameters from fitting
    fit = fits.open(fit_file)
    param = fit['PARAMETERS'].data

    # Covert wavelength to frequency since parameters are in frequency, and put exposure time into array
    freq = const.c.to('um/s').value / fitting_data['WAVELENGTH']
    exposure_time = fitting_data['EXPOSURE_MAP']

    # Find indices in frequency array that correspond to the peak frequency of the Guassian
    #   and its FWHM range. These indices will be used to determine the effective exposure time.
    hwhm = np.sqrt(2 * np.log(2)) * param[2, :, :]
    effect_exp_time = np.zeros(fitting_data['FLUX'].shape[1:3])
    for i in range(fitting_data['FLUX'].shape[1]):
        for j in range(fitting_data['FLUX'].shape[2]):
            if not np.isnan(param[1, i, j]):
                peak_ind = np.argmin(np.abs(freq - param[1, i, j]))
                effect_exp_time_ind = slice(np.argmin(np.abs(freq - freq[peak_ind] - hwhm[i, j])),
                                            np.argmin(np.abs(freq - freq[peak_ind] + hwhm[i, j]))+1)
                effect_exp_time[i, j] = np.nanmean(exposure_time[effect_exp_time_ind, i, j])

    # Construct primary header
    primary_hdu = fits.PrimaryHDU(header=primary_hdr)
    # Construct the image extensions and assign appropriate units and extension names
    image_hdr['BUNIT'] = ''
    image_hdr['EXTNAME'] = 'EFFECTIVE_EXPOSURE_TIME'
    exptime_hdu = fits.ImageHDU(effect_exp_time, header=image_hdr)
    # Create hdu list and save to output file
    new_hdul = fits.HDUList([primary_hdu, exptime_hdu])
    new_hdul.writeto(out_filename, overwrite=True)
