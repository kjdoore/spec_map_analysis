def gauss_func(x, a0, a1, a2, a3=None, a4=None, a5=None):
    """
    Defines a function that consists of a Gaussian with the optional
     addition of a polynomial up to degree 2.
    Parameters
    ----------
    x : 1-D array_like
        The independent variable data, of length M, where the function
        is to be evaluated.
    a0 : scalar
        The parameter that gives the height of the Gaussian in the
        function given by the equation:
        f = a0*exp(-((x-a1)/a2)^2/2) + a3 + a4*x + a5*x^2
    a1 : scalar
        The parameter that gives the location of the center of the
        Gaussian in the above equation.
    a2 : scalar
        The parameter that gives the sigma (width) of the Gaussian in
        the above equation.
    a3 : scalar, optional
        The parameter that gives the constant polynomial term in the
        above equation. If not specified, then no constant term is
        included in the function.
    a4 : scalar, optional
        The parameter that gives the linear polynomial term in the
        above equation. If not specified, then no linear term is
        included in the function.
    a5 : scalar, optional
        The parameter that gives the quadratic polynomial term in
        the above equation. If not specified, then no quadratic term
        is included in the function.
    Returns
    -------
    fx : 1-D array
        The dependent variable data, of length M, as determined from
        the above function for each value in x
    """

    import numpy as np

    # Make x a numpy array
    x = np.array(x)

    # Determine the number of terms to include in
    #   the function
    nterms = 3
    if a3 is not None:
        nterms = 4
    if a4 is not None:
        nterms = 5
    if a5 is not None:
        nterms = 6

    # Check to make sure the width is non-zero and positive.
    #   If it is not, then assume there is no Gaussian.
    if a2 > 0:
        z = (x - a1) / a2
        fx = a0 * np.exp(-z ** 2 / 2)
    else:
        fx = np.repeat(0, x.size)

    # Generate the resulting output
    if nterms == 4:
        fx = fx + a3
    elif nterms == 5:
        fx = fx + a3 + a4 * x
    elif nterms == 6:
        fx = fx + a3 + a4 * x + a5 * x ** 2

    return fx


def gaussfunc_jacobian(x, a0, a1, a2, a3=None, a4=None, a5=None):
    """
    Defines the Jacobian matrix of a Gaussian with the optional
      addition of a polynomial up to degree 2 with respect to the
      parameters.
    Parameters
    ----------
    x : 1-D array_like
        The independent variable data, of length M, where the Jacobian
        is to be evaluated.
    a0 : scalar
        The parameter that gives the height of the Gaussian in the
        function given by the equation:
        f = a0*exp(-((x-a1)/a2)^2/2) + a3 + a4*x + a5*x^2
    a1 : scalar
        The parameter that gives the location of the center of the
        Gaussian in the above equation.
    a2 : scalar
        The parameter that gives the sigma (width) of the Gaussian
        in the above equation.
    a3 : scalar, optional
        The parameter that gives the constant polynomial term in
        the above equation. If not specified, then the Jacobian
        does not include this parameter.
    a4 : scalar, optional
        The parameter that gives the linear polynomial term in
        the above equation. If not specified, then the Jacobian
        does not include this parameter.
    a5 : scalar, optional
        The parameter that gives the quadratic polynomial term in
        the above equation. If not specified, then the Jacobian
        does not include this parameter.
    Returns
    -------
    dfx : 2-D array
        The Jacobian matrix, an (M, k)-shaped array, of the above
        function for each value in x, where k is the number of
        function parameters specified
    """

    import numpy as np

    # Make x a numpy array
    x = np.array(x)

    # Determine the number of terms to include in
    #   the function
    nterms = 3
    if a3 is not None:
        nterms = 4
    if a4 is not None:
        nterms = 5
    if a5 is not None:
        nterms = 6

    # Check to make sure width is non-zero and positive.
    #   If it is not then assume there is no Gaussian.
    #   Compute partial derivatives with respect to each parameter.
    if a2 > 0:
        z = (x - a1) / a2
        d0 = np.exp(-z ** 2 / 2)
        d1 = a0 * z / a2 * d0
        d2 = d1 * z
    else:
        d0 = np.repeat(0, x.size)
        d1 = np.repeat(0, x.size)
        d2 = np.repeat(0, x.size)

    # Combine derivatives into Jacobian matrix
    if nterms == 3:
        dfx = np.column_stack((d0, d1, d2))
    elif nterms == 4:
        d3 = np.repeat(1.0, np.size(x))
        dfx = np.column_stack((d0, d1, d2, d3))
    elif nterms == 5:
        d3 = np.repeat(1.0, np.size(x))
        d4 = x
        dfx = np.column_stack((d0, d1, d2, d3, d4))
    elif nterms == 6:
        d3 = np.repeat(1.0, np.size(x))
        d4 = x
        d5 = x ** 2
        dfx = np.column_stack((d0, d1, d2, d3, d4, d5))

    return dfx


def gaussfunc_best_guess(x, y, nterms=3, min_bounds=None, max_bounds=None):
    """
    Generates and initial guess of each parameter when fitting
      a Gaussian with the optional addition of a polynomial up
      to degree 2 using non-linear least squares methods.
    Parameters
    ----------
    x : 1-D array_like
        The independent variable data, of length M
    y : 1-D array_like
        The dependent variable data, of length M, which a Gaussian
        function with the optional addition of a polynomial up
        to degree 2 is to be fit.
    nterms : integer, optional
        An integer from 3 to 6 specifying the number of terms
        to include in the Gaussian function given by:
        f = a0*exp(-((x-a1)/a2)^2/2) + a3 + a4*x + a5*x^2
        If only a value of 3 is specified, then only estimates for a
        Gaussian are returned (a0, a1, a2). If a value > 3 is specified
        then a polynomial of degree = nterms - 4 is added to the
        Gaussian up to a degree of 2 (a3, a4, a5).
        If a value is not specified, then a default value of 3 is used.
    min_bounds : 1-D array-like, optional
        An array of length nterms giving the minimum bounds for each
        parameter to be used in the non-linear least squares fitting.
        Guesses are restricted to be larger than these values.
    max_bounds : 1-D array-like, optional
        An array of length nterms giving the maximum bounds for each
        parameter to be used in the non-linear least squares fitting.
        Guesses are restricted to be smaller than these values.
    Returns
    -------
    parameter_guess : 1-D array
        An array of length nterms giving the initial guesses for each
        parameter to be used in the non-linear least squares fitting.
    """

    import numpy as np

    # Check to make sure inputs are of correct form
    # x and y must be 1-D arrays and same size
    x = np.array(x)
    y = np.array(y)
    if x.ndim != 1 or y.ndim != 1:
        raise ValueError("x and y must be one dimensional")
    if x.size != y.size:
        raise ValueError("x and y must be the same length")
    # nterms needs to be an integer between 3 and 6. So, check if integer and check range.
    nterms = int(nterms)
    if not 3 <= nterms <= 6:
        raise ValueError(
            ("nterms must be between 3 and 6; value given: {0}"
             ).format(nterms)
        )
    # min_bounds and max_bounds need to be nterms long if specified else give
    #   them infinite range
    if min_bounds is not None:
        if len(min_bounds) != nterms:
            raise ValueError("min_bounds must have nterms number of elements")
    else:
        min_bounds = np.repeat(-np.inf, nterms)
    if max_bounds is not None:
        if len(max_bounds) != nterms:
            raise ValueError("max_bounds must have nterms number of elements")
    else:
        max_bounds = np.repeat(np.inf, nterms)
    # min_bounds must be smaller or equal to max_bounds
    if any(min_bounds > max_bounds):
        raise ValueError("max_bounds must be larger than min_bounds")

    # For a Gaussian with a polynomial, subtract off a constant or line
    #   to get good initial estimate. Use a constant if only a constant
    #   term is used, and a line if 1st or 2nd order terms are used.
    #   No quadratic fit is used even with a 2nd order term, due to the
    #   Gaussian and quadratic potentially being highly correlated. If
    #   the polynomial coefficients are outside the parameter bounds,
    #   set them just inside the bounds.
    if nterms == 4:
        deg = 0
        coeff = np.polynomial.polynomial.polyfit(x, y, deg)
        if coeff[0] < min_bounds[3]:
            coeff[0] = min_bounds[3] + 1e-6
        if coeff[0] > max_bounds[3]:
            coeff[0] = max_bounds[3] - 1e-6
        gauss_res = y - np.polynomial.polynomial.polyval(x, coeff)
    elif nterms >= 5:
        deg = 1
        coeff = np.polynomial.polynomial.polyfit(x, y, deg)
        if coeff[0] < min_bounds[3]:
            coeff[0] = min_bounds[3] + 1e-6
        if coeff[0] > max_bounds[3]:
            coeff[0] = max_bounds[3] - 1e-6
        if coeff[1] < min_bounds[4]:
            coeff[1] = min_bounds[4] + 1e-6
        if coeff[1] > max_bounds[4]:
            coeff[1] = max_bounds[4] - 1e-6
        gauss_res = y - np.polynomial.polynomial.polyval(x, coeff)
    else:
        gauss_res = y

    # Find the min and max values of the data after removing constant
    #   or line (residual). Force the peak values to be within the
    #   min/max bounds of the peak offset (a1).
    gauss_res_limit = np.where(np.logical_or(x < min_bounds[1], x > max_bounds[1]), np.nan, gauss_res)
    ymax = np.nanmax(gauss_res_limit)
    ymin = np.nanmin(gauss_res_limit)

    # Find the largest deviation from the residual and its location
    #   in the independent variable space.
    if np.abs(ymax) > np.abs(ymin):
        peak_loc = np.nanargmax(gauss_res_limit)
    else:
        peak_loc = np.nanargmin(gauss_res_limit)
    # If the largest deviation is outside the parameter bounds,
    #   then set to the other deviation
    if min_bounds[0] > ymin and max_bounds[0] < ymax:
        raise ValueError("Peak height bounds are to restrictive to find a peak")
    elif min_bounds[0] > ymin:
        peak_loc = np.nanargmax(gauss_res_limit)
    elif max_bounds[0] < ymax:
        peak_loc = np.nanargmin(gauss_res_limit)

    # Do not allow the peak location to be on the edge of the
    #   independent variable array. If it is, move it in one
    #   index.
    if peak_loc == 0:
        peak_loc = 1
    elif peak_loc == (y.size - 1):
        peak_loc = y.size - 2

    # Set the peak guess as the peak value of the residual. If
    #   the estimated peak height is larger than the maximum bound,
    #   set to just below maximum bound. If it is smaller than the
    #   minimum bound, set just above.
    peak_height = gauss_res[peak_loc]
    if peak_height > max_bounds[0]:
        peak_height = max_bounds[0] - 1e-6
    if peak_height < min_bounds[0]:
        peak_height = min_bounds[0] + 1e-6

    # The expected height of the Gaussian at the one sigma
    #   deviation from the peak/trough.
    sig_height = gauss_res[peak_loc] / np.exp(1)

    # Estimate the location of the one sigma height of
    #   the Gaussian by finding the independent variable index
    #   where the dependent variable is just less than the
    #   one sigma height.
    i = 0
    while ((peak_loc + i + 1) < y.size) and ((peak_loc - i) > 0) and \
            (np.abs(gauss_res[peak_loc + i]) > np.abs(sig_height)) and \
            (np.abs(gauss_res[peak_loc - i]) > np.abs(sig_height)):
        i += 1
    # The sigma is the difference between the peak location and
    #  the one sigma height location.
    sigma = np.abs(x[peak_loc] - x[peak_loc + i])
    # If the estimated sigma is larger than the maximum bound,
    #   set to just below maximum bound. If it is smaller than the
    #   minimum bound, set just above.
    if sigma > max_bounds[2]:
        sigma = max_bounds[2] - 1e-6
    if sigma < min_bounds[2]:
        sigma = min_bounds[2] + 1e-6

    # Generate the initial best guess for the Gaussian parameters.
    #   The peak location is the location of the peak.
    parameter_guess = np.array([peak_height, x[peak_loc], sigma])

    # If a polynomial is added to the Gaussian include those estimates,
    #   which are the coefficients from the constant/line subtraction.
    if nterms > 3:
        parameter_guess = np.append(parameter_guess, coeff[0])

    if nterms > 4:
        parameter_guess = np.append(parameter_guess, coeff[1])

    # Assume the quadratic term is zero to prevent correlation with
    #   Gaussian.
    if nterms > 5:
        parameter_guess = np.append(parameter_guess, 0)

    return parameter_guess


def gaussfunc_bounds(x, y, nterms=3):
    """
    Generates bounds based on the data for each parameter when
      fitting a Gaussian with the optional addition of a polynomial
      up to degree 2 using non-linear least squares methods.
    Parameters
    ----------
    x : 1-D array_like
        The independent variable data, of length M
    y : 1-D array_like
        The dependent variable data, of length M, which a Gaussian
        function with the optional addition of a polynomial up
        to degree 2 is to be fit
    nterms : integer, optional
        An integer from 3 to 6 specifying the number of terms
        to include in the Gaussian function given by:
        f = a0*exp(-((x-a1)/a2)^2/2) + a3 + a4*x + a5*x^2
        If only a value of 3 is specified then only bounds for a
        Gaussian are returned (a0, a1, a2). If a value > 3 is specified
        then bounds for the polynomial terms of degree = nterms - 1
        up to a degree of 2 (a3, a4, a5) are included.
        If a value is not specified, then a default value of 3 is used.
    Returns
    -------
    min_bounds : 1-D array
        An array of length nterms giving the minimum bounds for each
        parameter to be used in the non-linear least squares fitting.
    max_bounds : 1-D array
        An array of length nterms giving the maximum bounds for each
        parameter to be used in the non-linear least squares fitting.
    """

    import numpy as np

    # Check to make sure inputs are of correct form
    # x and y must be 1-D arrays and same size
    x = np.array(x)
    y = np.array(y)
    if x.ndim != 1 or y.ndim != 1:
        raise ValueError("x and y must be one dimensional")
    if x.size != y.size:
        raise ValueError("x and y must be the same length")
    # nterms needs to be an integer between 3 and 6. So, make integer and check range.
    nterms = int(nterms)
    if not 3 <= nterms <= 6:
        raise ValueError(
            ("nterms must be between 3 and 6; value given: {0}"
             ).format(nterms)
        )

    # The height bound of the Gaussian cannot be larger than 1.5x the
    #   min/max y value. This allows for the Gaussian to be slightly
    #   larger than the min/max value. The center of the
    #   Gaussian cannot be outside the x value range, and the width
    #   cannot be negative or larger than the range of x.
    min_bounds = np.array([np.min(y) * 1.5, np.min(x), 0])
    max_bounds = np.array([np.max(y) * 1.5, np.max(x), np.ptp(x)])

    # The constant polynomial term cannot be larger/smaller than the max/min
    #   y value.
    if nterms >= 4:
        min_bounds = np.append(min_bounds, np.min(y))
        max_bounds = np.append(max_bounds, np.max(y))

    # The linear and quadratic terms have no bounds.
    if nterms >= 5:
        min_bounds = np.append(min_bounds, -np.inf)
        max_bounds = np.append(max_bounds, np.inf)
    if nterms == 6:
        min_bounds = np.append(min_bounds, -np.inf)
        max_bounds = np.append(max_bounds, np.inf)

    return min_bounds, max_bounds
