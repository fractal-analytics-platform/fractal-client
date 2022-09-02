"""
Copyright 2022 (C) Friedrich Miescher Institute for Biomedical Research and
University of Zurich

Original authors:
Tommaso Comparin <tommaso.comparin@exact-lab.it>
Marco Franzon <marco.franzon@exact-lab.it>

This file is part of Fractal and was originally developed by eXact lab S.r.l.
<exact-lab.it> under contract with Liberali Lab from the Friedrich Miescher
Institute for Biomedical Research and Pelkmans Lab from the University of
Zurich.
"""
import shutil

import dask.array as da


def to_zarr_custom(newzarrurl=None, array=None, component="", overwrite=False):

    """

    Custom workaround for dask issue
    (https://github.com/dask/dask/issues/5942), where a dask array loaded with
    from_zarr cannot be written with to_zarr(..., overwrite=True).

    :param newzarrurl: output zarr file
    :type newzarrurl: str
    :param array: dask array to be stored
    :type array: dask array
    :param component: target subfolder of the zarr file (optional, default "")
    :type component: str
    :param overwrite: overwrite existing data (optional, default False)
    :type overwrite: boolean

    """

    # Check required arguments
    if newzarrurl is None:
        raise Exception("ERROR: Missing newzarrurl arg in to_zarr_custom")
    if array is None:
        raise Exception("ERROR: Missing array arg in to_zarr_custom")

    # Sanitize arguments
    if component.endswith("/"):
        component = component[:-1]
    if not newzarrurl.endswith("/"):
        newzarrurl += "/"

    if overwrite:
        tmp_suffix = "_TEMPORARY_ZARR_ARRAY"
        array.to_zarr(
            newzarrurl,
            component=component + tmp_suffix,
            dimension_separator="/",
            compute=True,
        )
        shutil.rmtree(newzarrurl + component)
        shutil.move(
            newzarrurl + component + tmp_suffix, newzarrurl + component
        )
        output = da.from_zarr(newzarrurl, component=component)
    else:
        output = array.to_zarr(
            newzarrurl,
            component=component,
            dimension_separator="/",
            compute=True,
            return_stored=True,
        )

    return output
