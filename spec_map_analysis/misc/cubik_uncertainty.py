def cubik_uncertainty(out_filename, pipeline_filename, cubik_filename, err_extension, **kwargs):
    """
    A function that reads a level 4 FIFI-LS pipeline FITS file
      and an output file from Cubik and replaces the pipeline
      errors with the Cubik errors. It also adds a HISTORY
      keyword stating this change.
      (Cubik source: https://github.com/darioflute/fifipy)
    Parameters
    ----------
    out_filename : string
        A string containing the updated file to write
    pipeline_filename : string
        A string containing the pipeline file to be opened
    cubik_filename : string
        A string containing the Cubik file to be opened
    err_extension : two element list of strings or integers
        A two element list of strings or integers containing
        the name or number of the error extensions for the
        pipeline and cubik files, respectively.
    kwargs
        Keyword arguments that get passed to the FITS
        writeto() function
    """

    from astropy.io import fits

    # Read in the pipeline and cubik files
    pipeline_hdu = fits.open(pipeline_filename)
    cubik_hdu = fits.open(cubik_filename)

    # Add a HISTORY keyword to the error extension header stating the replacement of errors
    cubik_hdu[err_extension[1]].header['HISTORY'] = 'Error estimated by Cubik (D. Fadda)'
    # Replace the pipeline error cube and header with the cubik error cube and header
    pipeline_hdu[err_extension[0]] = cubik_hdu[err_extension[1]]
    # Add a HISTORY keyword to the primary header stating the replacement of errors
    pipeline_hdu[0].header['HISTORY'] = 'Replaced uncertaities with Cubik estimate'
    # Write the file
    pipeline_hdu.writeto(out_filename, **kwargs)
