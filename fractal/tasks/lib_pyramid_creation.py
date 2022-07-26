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
import dask.array as da
import numpy as np

from fractal.tasks.lib_to_zarr_custom import to_zarr_custom


def write_pyramid(
    data,
    overwrite=False,
    newzarrurl=None,
    coarsening_xy=2,
    num_levels=2,
    chunk_size_x=None,
    chunk_size_y=None,
    aggregation_function=None,
):

    """
    Take a four-dimensional array and build a pyramid of coarsened levels

    :param data_czyx: input data
    :type data_czyx: dask array
    :param coarsening_xy: coarsening factor along X and Y
    :type coarsening_xy: int
    :param num_levels: number of levels in the zarr pyramid
    :type num_levels: int
    :param chunk_size_x: chunk size along X
    :type chunk_size_x: int
    :param chunk_size_y: chunk size along Y
    :type chunk_size_y: int
    :param aggregation_function: FIXME
    :type aggregation_function: FIXME
    """

    # Check the number of axes and identify YX dimensions
    ndims = len(data.shape)
    if ndims not in [2, 3, 4]:
        raise Exception(
            "ERROR: data has shape {data.shape}, ndims not in [2,3,4]"
        )
    y_axis = ndims - 2
    x_axis = ndims - 1

    # Set rechunking options, if needed
    if chunk_size_x is None or chunk_size_y is None:
        apply_rechunking = False
    else:
        apply_rechunking = True
        chunking = {y_axis: chunk_size_y, x_axis: chunk_size_x}

    # Set aggregation_function
    if aggregation_function is None:
        aggregation_function = np.mean

    # Write highest-resolution level
    level0 = to_zarr_custom(
        newzarrurl=newzarrurl, array=data, component="0", overwrite=overwrite
    )
    if apply_rechunking:
        previous_level = level0.rechunk(chunking)
    else:
        previous_level = level0

    # Compute and write lower-resolution levels
    for ind_level in range(1, num_levels):

        # Verify that coarsening is doable
        if min(previous_level.shape[-2:]) < coarsening_xy:
            raise Exception(
                f"ERROR: at {ind_level}-th level, "
                f"coarsening_xy={coarsening_xy} "
                f"but previous level has shape {previous_level.shape}"
            )

        # Apply coarsening
        newlevel = da.coarsen(
            aggregation_function,
            previous_level,
            {y_axis: coarsening_xy, x_axis: coarsening_xy},
            trim_excess=True,
        ).astype(data.dtype)

        # Apply rechunking
        if apply_rechunking:
            newlevel_rechunked = newlevel.rechunk(chunking)
        else:
            newlevel_rechunked = newlevel

        # Write zarr and store output (useful to construct next level)
        previous_level = to_zarr_custom(
            newzarrurl=newzarrurl,
            array=newlevel_rechunked,
            component=f"{ind_level}",
            overwrite=overwrite,
        )
