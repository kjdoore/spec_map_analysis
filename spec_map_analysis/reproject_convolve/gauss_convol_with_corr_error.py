def gauss_convol_with_corr_error(data, unc, sigma_smooth, sigma_beam, **kwargs):
    """
    Convolves the data with a Gaussian of the given sigma. The
      uncertainties are propagated assuming the uncertainties are
      correlated like a Gaussian convolution with the given beam size.
      Method is presented in Klein, R. 2021 RNAAS 5 39
      (DOI: 10.3847/2515-5172/abe8df)
    Parameters
    ----------
    data : array-like
        A one or two dimensional array giving the data that is to be
        convolved with the Gaussian
    unc : array-like
        A one or two dimensional array giving the uncertainties on the
        data that are to be propagated
    sigma_smooth : scalar
        A scalar giving the standard deviation value for the Gaussian
        convolution in pixels.
    sigma_beam : scalar
        A scalar giving the standard deviation value for the estimated
        correlation length (due to a beam size) in pixels.
    kwargs
        Keyword arguments passed to astropy.convolution.Gaussian2DKernel
    Returns
    -------
    convol_data : array
        A one or two dimensional array giving the convolved data
    new_unc : array
        A one or two dimensional array giving the propagated uncertainties
        of the convolved data
    """

    import numpy as np
    from astropy.convolution import convolve
    from astropy.convolution import Gaussian2DKernel
    from astropy.convolution import Gaussian1DKernel

    n_dim = data.ndim

    # Generate a Gaussian kernel depending on the data dimensions
    if n_dim == 1:
        gauss_kernel = Gaussian1DKernel(sigma_smooth, **kwargs)
    elif n_dim == 2:
        gauss_kernel = Gaussian2DKernel(sigma_smooth, **kwargs)
    else:
        raise ValueError("Input data must be one or two dimensional")

    # Convolve the data preserving the NaNs in the original data
    convol_data = convolve(data, gauss_kernel, preserve_nan=True, nan_treatment='fill')

    # Convolve the uncertainties with the non-normalized square of the kernel
    varianz = convolve(unc ** 2, gauss_kernel.array ** 2, normalize_kernel=False, nan_treatment='fill',
                       preserve_nan=True)

    # Generate new uncertainties to account for the correlation between data points.
    new_unc = (2 * np.sqrt(np.pi) * sigma_smooth * sigma_beam /
               np.sqrt(sigma_smooth**2 + sigma_beam**2)) ** (n_dim / 2.) * np.sqrt(varianz)

    return convol_data, new_unc
