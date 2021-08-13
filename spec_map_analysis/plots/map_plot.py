def map_plot(out_filename, input_file0, extension0, cmap='viridis', ralim=None, declim=None,
             figure_height_scale=1, figure_width_scale=1, **kwargs):
    """
    A versatile image plotting function that can plot any number of
      input images, with the only requirement being that the images
      are of the same or nearby objects (i.e., They must be within
      a similar RA and Dec range. Large differences in RA and Dec
      compared to the pixel scale will result in tiny indistinguishable
      images begin plotted.). See kwargs for details on plotting
      multiple images in the same figure, along with other functionality.
    Parameters
    ----------
    out_filename : string
        A string containing the file that the figure is to be saved.
    input_file0 : string
        A string containing the FITS file to be read.
    extension0 : string or int
        A string or integer containing the name or index of
        the extension of input_file0 that is to be read.
    cmap : string, optional
        A string containing the colormap that the images are
        to be plotted.
    ralim : list of strings, optional
        A two element list containing the upper and lower Right Ascension
        limits as strings to be used when plotting. Values must be
        clearly discernible by the astropy.coordinates.Angle function.
        For examples see:
        https://docs.astropy.org/en/stable/api/astropy.coordinates.Angle.html#astropy.coordinates.Angle
        If not specified, then the RA limits are automatically determined
        from the input_file0 image.
    declim : list of strings
        A two element list containing the upper and lower Declination
        limits as strings to be used when plotting. Values must be
        clearly discernible by the astropy.coordinates.Angle function.
        For examples see:
        https://docs.astropy.org/en/stable/api/astropy.coordinates.Angle.html#astropy.coordinates.Angle
        If not specified, then the Dec limits are automatically determined
        from the input_file0 image.
    figure_height_scale : scalar
        A scalar containing a value to scale the height of the
        output plot.
    figure_width_scale : scalar
        A scalar containing a value to scale the width of the
        output plot.
    kwargs
        In order to allow for an unlimited number of images to be
        plotted and for each to have the same misc functionality,
        any additional images beyond the single required image will
        be passed to kwargs, along with all functionality exclusive
        to a single image. A detailed description of all other
        optional accepted inputs are specified below.

        Note that the X at the end of all the keywords refers to
        a variable integer value, where X is the Xth image in python
        indexing. These keywords reflect functionality that applies
        to the Xth image. For example, to plot two additional images
        besides the required image, X would be 1 and 2 for the second and
        third image, respectively. If X values are not sequential
        (i.e., X is 1 and 3 rather than 1 and 2), all inputs after
        the skipped value will be ignored. This means that X values
        must be sequential with no skipping of integer values.

        input_fileX: string
            A string containing the Xth FITS file to be read.
        extensionX: string or integer
            A string or integer containing the name or index of
            the extension of input_fileX that is to be read.
        operatorX : list of strings
            A list of string(s) containing the operation that
            is to be performed to generate the plotted Xth image.

            This keyword allows for mathematical operations to be
            performed on the Xth input image. The code evoked when
            this keyword is given makes use of the built-in python
            function, eval(). This function takes a string and
            performs the python statement in the string.
            For example, eval('2+2') would output a value of 4.

            If using this keyword, the Xth image is read-in as the
            variable 'a'. Then, any python code contained in
            operatorX referencing 'a' will be performed on the image.
            So for example, if the logarithm of the input image was to
            be plotted, this keyword could be set to ['np.log10(a)'].
            This would cause the Xth plotted image to be the logarithm of the
            Xth input image. Another example, if the inverse of the input
            image was to be plotted, this keyword could be set to ['1/a'].
            This keyword allows for operations to be performed in situ
            rather than having to be performed outside of this function
            and saved to a new image extension in a FITS file.

            Further, one may want to plot images that require more complex
            operations combining multiple images without having to save
            the operation to a new image extension. This keyword allows
            for multiple additional images to be read-in, combined with
            the Xth input image, and the result plotted. These additional
            images and their extensions are read-in using the operator_fileX
            and operator_extensionX keywords (see below).

            Each additional image in the operator_fileX and operator_extensionX
            lists are consecutively read-in as the variable 'b'. Then, any
            python code contained in operatorX list referencing 'a' and 'b'
            will be consecutively performed and the result with be placed
            into the variable 'a'. For a simple example, if the ratio of
            two images was to be plotted, this keyword could be set to
            ['b/a'] or ['a/b'], depending on which image is to be the
            numerator/denominator. For a more complex example, if four
            images were to be combined with addition and division as
            w/x + y/z, where each image is a unique letter, the
            keyword could be given as a list that preforms the required
            operations using only two images at a time. To do this, rearrange
            the above expression to: (w*z/x + y)/z. Then operatorX can be
            set to ['a*b', 'a/b', 'a+b', 'a/b'], assuming w is the input
            image. Here, w is the original 'a' variable. It is then multiplied
            by z, which is read in as 'b', and the new array is updated as 'a'.
            Then this new 'a' is divided by x, which is again read in
            as 'b', and the new array is again updated as 'a'. This process
            is repeated for each index in the operatorX keyword. In other word,
            operations need to be structured in a fashion that allows for
            only two or less variables to be combined at a time, where
            'a' is the current image, and 'b' is the newly read in image
            from the operator_fileX and operator_extensionX list.

            See the examples for a complete demonstration on how to use
            this keyword.

            Notes: 1) Highly complex operations are allowed. However, each
                       operation can only contain at most two images,
                       must be wholly contained in a single string, not
                       reference any other variables besides 'a' and 'b',
                       and only use built-in python functions and numpy
                       functions.
                   2) If an operation is performed only on 'a', it
                       does not require a corresponding index in the
                       operator_fileX and operator_extensionX lists if
                       it is the last operation or only operation performed.
        operator_fileX : list of strings
            A list of string(s) containing the FITS file(s) to be
            read that will have an operation performed on input_fileX.
            See operatorX kwarg input for details.
        operator_extensionX : list of strings or integers
            A list of string(s) or integer(s) containing the name or
            index of the extension(s) of operator_fileX(s) that is to
            be read. See operatorX kwarg input for details.
        rangeX : two element array-like
            A two element array-like class that contains the color
            bar value range for the Xth image
        outlineX : integer
            An integer containing the image number whose outline is
            to be plotted on the Xth image.

            This keyword allows for an outline of an input image to
            be placed on the Xth image. This is convenient when images
            of the same object are not the same size and an outline
            of the smaller image is wanted on the larger image for
            easier comparison.

            Example: outline2=1
                    Means an outline of image 1 (2nd image in the
                    python notation used) is placed on image 2
                    (3rd image)
            Example: outline1=1
                    Means an outline of image 1 (2nd image in the
                    python notation used) is placed on itself.
        outline_edge_colorX : string
            A string containing the color of the outline on the
            Xth image (default: red)
        apertureX : list of two element lists of strings
            A list of two element lists of strings containing the
            Right Ascension and Declination location to place a
            circular aperture on the Xth image. Each set of two
            element list must have a corresponding radii given
            in aperture_radiusX.
            Values must be clearly discernible by the
            astropy.coordinates.Angle function.
            For examples see:
            https://docs.astropy.org/en/stable/api/astropy.coordinates.Angle.html#astropy.coordinates.Angle
        aperture_colorX : string
            A string containing the color of all apertures on the
            Xth image (default: red)
        aperture_radiusX : list of scalars
            A list of scalars containing the radii of the apertures
            in arcseconds on the Xth image
        titleX : string
            A string containing the title for the Xth image
        cbar_labelX : string
            A string containing the color bar label for the
            Xth image
        negativeX : string
            A string containing the color map that negative values
            in the image should be set.
        contourX : list of strings
            A list of string(s) containing the FITS file(s) to be
            read for contours to placed on the Xth image
        contour_extensionX : list of strings
            A list of string(s) or integer(s) containing the name or index
            of the extension(s) of each contourX file to be read.
        contour_colorX : list of strings
            A list of string(s) containing the colors for each set of
            contours on the Xth image. (default: white)
        contour_alphaX : list of scalars
            A list of scalars from 0 to 1 giving the opacity of each
            set of contours on the Xth image. (default: 1)
        contour_levelsX : list of lists of scalars
            A list of list of scalars containing the contour levels for
            each set of contours on the Xth image. (default: 3 contour
            values are linearly generated between the min and max data
            range (np.linspace(min,max,5)[1:4]))
    Notes
    -----
    - When plotting, additional images are plotted row first, left
        to right, then column, top to bottom.
    Examples
    -------
    1) This example shows how to plot flux and uncertainty images
       with titles and color bar labels from the same file with extension
       names 'FLUX' and 'UNCERTAINTY'.

        from plots.map_plot import map_plot
        file = <path_to_file>+'file.fits'
        ofile = <path_to_save_output_image>+'plot.pdf'
        map_plot(ofile, file, 'FLUX', title0='Flux', cbar_label0='Jy',
                 input_file1=file, extension1='UNCERTAINTY', title1='Flux Uncertainty',
                 cbar_label1='Jy')

    2) This example is the same as example (1) except it plots the flux and
       the normalized signal-to-noise ratio from the same file using
       the operator keywords.

        from plots.map_plot import map_plot
        file = <path_to_file>+'file.fits'
        ofile = <path_to_save_output_image>+'plot.pdf'
        map_plot(ofile, file, 'FLUX', title0='Flux', cbar_label0='Jy',
                 input_file1=file, extension1='UNCERTAINTY', title1='Normalized Signal-to-Noise',
                 cbar_label1='Normalized SNR', operator1=['b/a', 'a/np.nanmax(a)'], operator_file1=[file],
                 operator_extension1=['FLUX'])
        b/a * np.nanmax(b/a)
        # Note: operator_file1 and operator_extension1 do not need a second index due to operator1[1]
        #   being the last operation and not acting on another image.
    """

    from astropy.wcs import WCS
    import matplotlib.pyplot as plt
    from astropy.io import fits
    import numpy as np
    from astropy.coordinates import Angle
    from astropy.visualization.wcsaxes import Quadrangle
    from astropy.visualization.wcsaxes import SphericalCircle
    from astropy import units as u

    # Read in the required input FITS file extension and place it in a HDU dictionary
    hdu_dict = {'0': (fits.open(input_file0))[extension0]}

    # Read in any other input FITS file extensions and place them in the HDU dictionary
    # Loops though all misc function inputs and checks for other input files
    hdu_number = 1
    for key in kwargs:
        if 'input_file'+str(hdu_number) in key:
            hdu_dict[str(hdu_number)] = (fits.open(kwargs[key]))[kwargs['extension'+str(hdu_number)]]
            hdu_number += 1

    # Create a WCS object for each input file extension and place them into a dictionary
    wcs = {}
    for wcs_number, key in enumerate(hdu_dict):
        wcs[str(wcs_number)] = WCS(hdu_dict[key].header)

    # If the RA and DEC limits are set, convert them to pixel units for plot axis limits
    #   If not specified, set the limits to be the footprint of the first input image.
    if ralim is not None and declim is not None:
        radeg = Angle(ralim, unit='degree')
        decdeg = Angle(declim, unit='degree')
        xlim, ylim = wcs['0'].all_world2pix(radeg, decdeg, 0)
    else:
        footprint = wcs['0'].calc_footprint()
        radeg = Angle([np.max(footprint[:, 0]), np.min(footprint[:, 0])], unit='degree')
        decdeg = Angle([np.min(footprint[:, 1]), np.max(footprint[:, 1])], unit='degree')
        xlim, ylim = wcs['0'].all_world2pix(radeg, decdeg, 0)

    # Based off of the number of input files, determine how many rows and columns are needed
    #   for subplot. The subplot grid shape will always be square or rectangular with
    #   ncols >= nrows
    num_plots = len(hdu_dict)
    nrows = int(np.floor(np.sqrt(num_plots)))
    if nrows != 1:
        ncols = int(np.ceil(num_plots/nrows))
    else:
        ncols = num_plots

    # Generate the figure with adjustable width and height
    fig = plt.figure(figsize=(ncols * figure_width_scale, nrows * figure_height_scale))

    # Loop through each input file and plot using the required file's WCS as the axis reprojection
    for index in range(0, num_plots):
        ax = fig.add_subplot(nrows, ncols, index + 1, projection=wcs['0'])

        # Generate the data to be plotted. If an operation is to be performed on the data,
        #   preform the operation. If extra data is needed to preform the operation, read in
        #   the extension of the file and preform the operation for each extra data file.
        a = hdu_dict[str(index)].data
        if 'operator' + str(index) in kwargs:
            for i, operation in enumerate(kwargs['operator' + str(index)]):
                if 'operator_file' + str(index) in kwargs and 'operator_extension' + str(index) in kwargs:
                    if len(kwargs['operator_file' + str(index)])-1 >= i:
                        b = (fits.open(kwargs['operator_file' + str(index)][i]))[kwargs['operator_extension' +
                                                                                        str(index)][i]].data
                a = eval(operation)
        data = a

        if 'negative'+str(index) in kwargs:
            # Copy the data and replace positive values with NaN to prevent covering positive data in first plot
            #   and replace negative data with NaNs in positive copy
            negative_data = data + 0
            data[data < 0] = np.nan
            negative_data[data > 0] = np.nan

        # Plot the image and limit the color range to the optional input values. If not the required
        #   image, transform the WCS of the image axis to the required image's WCS axis
        if index == 0:
            if 'range'+str(index) in kwargs:
                im = ax.imshow(data, origin='lower', vmin=(kwargs['range' + str(index)])[0],
                               vmax=(kwargs['range' + str(index)])[1], cmap=cmap)
            else:
                im = ax.imshow(data, origin='lower', cmap=cmap)
            # If negative keyword is specified, plot the negative data as a different color map over current data.
            if 'negative'+str(index) in kwargs:
                im2 = ax.imshow(negative_data, origin='lower', cmap=kwargs['negative'+str(index)], alpha=1.0)
        else:
            if 'range' + str(index) in kwargs:
                im = ax.imshow(data, origin='lower', vmin=(kwargs['range' + str(index)])[0], cmap=cmap,
                               vmax=(kwargs['range' + str(index)])[1], transform=ax.get_transform(wcs[str(index)]))
            else:
                im = ax.imshow(data, origin='lower', transform=ax.get_transform(wcs[str(index)]), cmap=cmap)
            if 'negative'+str(index) in kwargs:
                im2 = ax.imshow(negative_data, origin='lower', cmap=kwargs['negative'+str(index)], alpha=1.0,
                                transform=ax.get_transform(wcs[str(index)]))

        # If an outline is specified, plot the specified outline over the current image
        if 'outline'+str(index) in kwargs:
            # Gets the corners of the specified image in RA and DEC
            foot_print = wcs[str(kwargs['outline' + str(index)])].calc_footprint()
            # Determines the span of the outline
            foot_print_range = np.ptp(foot_print, axis=0)
            # Set the color of the outline if specified
            if 'outline_edge_color'+str(index) in kwargs:
                edge_color = kwargs['outline_edge_color'+str(index)]
            else:
                edge_color = 'red'
            # Generates the outline using the corners and span in units of degrees, assuming the corners
            #   are RA and DEC in fk5
            outline = Quadrangle(np.min(foot_print, axis=0) * u.deg, foot_print_range[0] * u.deg,
                                 foot_print_range[1] * u.deg, edgecolor=edge_color, facecolor='none',
                                 transform=ax.get_transform('fk5'))
            ax.add_patch(outline)

        # If an aperture is specified, plot the aperture over the current image
        if 'aperture'+str(index) in kwargs:
            for i, pos in enumerate(kwargs['aperture'+str(index)]):
                # Gets the center of the aperture image in RA and DEC in degrees
                position = [Angle(pos[0], unit='degree'), Angle(pos[1], unit='degree')]
                # Set the color of the aperture if specified
                if 'aperture_color'+str(index) in kwargs:
                    edge_color = kwargs['aperture_color'+str(index)]
                else:
                    edge_color = 'red'
                # Generates the aperture using the center position and radius in units of arcseconds,
                #   assuming the center RA and DEC in fk5
                aper = SphericalCircle(position, (kwargs['aperture_radius'+str(index)])[i] * u.arcsec,
                                       edgecolor=edge_color, facecolor='none',
                                       transform=ax.get_transform('fk5'))
                ax.add_patch(aper)

        # If a contour is specified, plot the contour over the current image
        if 'contour'+str(index) in kwargs:
            for i, cont in enumerate(kwargs['contour'+str(index)]):
                contour_hdu = (fits.open(cont))[kwargs['contour_extension'+str(index)][i]]
                # Set the color of the contour if specified
                if 'contour_color'+str(index) in kwargs:
                    contour_color = kwargs['contour_color'+str(index)][i]
                else:
                    contour_color = 'white'
                # Set the opacity of the contour
                if 'contour_alpha'+str(index) in kwargs:
                    contour_alpha = kwargs['contour_alpha'+str(index)][i]
                else:
                    contour_alpha = 1
                # Generate 3 contour levels if not specified
                if 'contour_levels'+str(index) in kwargs:
                    contour_levels = kwargs['contour_levels'+str(index)][i]
                else:
                    contour_levels = np.linspace(np.nanmin(contour_hdu.data), np.nanmax(contour_hdu.data), 5)[1:4]
                # Generate the contours
                ax.contour(contour_hdu.data, levels=contour_levels, colors=contour_color, alpha=contour_alpha,
                           transform=ax.get_transform(WCS(contour_hdu.header)))

        # Set the axis limits if specified and label them
        ax.set(xlim=xlim, ylim=ylim, xlabel='RA', ylabel='Dec')
        ra = ax.coords[0]
        ra.set_major_formatter('hh:mm:ss')
        dec = ax.coords[1]
        dec.set_major_formatter('dd:mm:ss')

        # If specified, give the plot a title
        if 'title'+str(index) in kwargs:
            ax.set_title(kwargs['title' + str(index)])

        # Add a colorbar to the plot and label it if specified
        cbar = fig.colorbar(im, ax=ax)
        if 'cbar_label' + str(index) in kwargs:
            cbar.set_label(kwargs['cbar_label' + str(index)])

    # Save the resulting figure
    fig.savefig(out_filename, dpi=500)
    plt.close()
