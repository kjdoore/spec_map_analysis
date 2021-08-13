def log_prob(param, x, y, uncertainty, lower_bounds=None, upper_bounds=None):
    """
    Function to compute the log-likelihood probability for input
      into emcee. Calculates the non-normalized probability of the
      quality of fit of the a Gaussian model to the data using the
      chisqr statistic.
    Parameters
    ----------
    param : 1-D array-like
        An array of length 3 to 6 containing the parameters for
        the Gaussian function given by:
        f = a0*exp(-((x-a1)/a2)^2/2) + a3 + a4*x + a5*x^2
    x : 1-D array-like
        An array containing the independent variable values
    y : 1-D array-like
        An array containing the measured dependent variable values
        at each x value
    uncertainty : 1-D array-like
        An array containing the uncertainties on the measured
        dependent variable values at each x value
    lower_bounds : 1-D array-like, optional
        An array with the same length as param giving the lower bounds
        beyond which the probability is 0
    upper_bounds : 1-D array-like, optional
        An array with the same length as param giving the upper bounds
        beyond which the probability is 0
    Returns
    -------
    log_prob : scalar
        The log-likelihood probability of the model to the data using
        the chisqr statistic
    """

    from spec_map_analysis.spectra_fitting import gauss_func as gf
    import numpy as np

    # Check to make sure inputs are of correct form
    # param must have 3 to 6 elements
    param = np.array(param)
    if not 3 <= param.size <= 6:
        raise ValueError(
            ("param must have between 3 and 6 elements; currently has {0} elements"
             ).format(param.size)
        )
    # x, y, and uncertainty must be 1-D arrays and same size
    x = np.array(x)
    y = np.array(y)
    uncertainty = np.array(uncertainty)
    if x.ndim != 1 or y.ndim != 1 or uncertainty.ndim != 1:
        raise ValueError("x, y, and/or uncertainty must be one dimensional")
    if x.size != y.size or x.size != uncertainty.size or y.size != uncertainty.size:
        raise ValueError("x, y, and uncertainty must be the same length")

    # Check if bounds are specified, if not set to +/- infinity
    if lower_bounds is None:
        lower_bounds = -np.inf
    if upper_bounds is None:
        upper_bounds = np.inf

    # If any parameters is outside of the bounds return a log-prob
    #   of -infinity, which is a prob of 0.
    if any(param > upper_bounds) or any(param < lower_bounds):
        return -np.inf

    # Calculate model data and chisqr
    mod_y = gf.gauss_func(x, *param)
    chisqr = np.nansum(((mod_y - y) / uncertainty) ** 2)

    return -0.5 * chisqr
