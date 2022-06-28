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


def create_pyramid(
    data_czyx,
    coarsening_z=1,
    coarsening_xy=1,
    num_levels=1,
    chunk_size_x=None,
    chunk_size_y=None,
    num_channels=None,
):

    """
    Take a four-dimensional array and build a pyramid of coarsened levels

    :param data_czyx: input data
    :type data_czyx: dask array
    :param coarsening_z: coarsening factor along Z
    :type coarsening_z: int
    :param coarsening_xy: coarsening factor along X and Y
    :type coarsening_xy: int
    :param num_levels: number of levels in the zarr pyramid
    :type num_levels: int
    :param chunk_size_x: chunk size along X
    :type chunk_size_x: int
    :param chunk_size_y: chunk size along Y
    :type chunk_size_y: int
    :param num_channels: number of channels
    :type num_channels: int
    """

    # Check that input has the right shape
    if len(data_czyx.shape) != 4:
        raise Exception("Error in create_pyramid: data_czyx has wrong shape")

    # Check rechunking options
    apply_rechunking = True
    if chunk_size_x is None or chunk_size_y is None:
        apply_rechunking = False

    # Coarsen globally along Z direction
    if coarsening_z > 1:
        data_czyx = da.coarsen(
            np.min, data_czyx, {1: coarsening_z}, trim_excess=True
        )

    # Create pyramid of XY-coarser levels
    pyramid = []
    for level in range(num_levels):
        list_data_channels = []
        for ind_chl in range(num_channels):

            # Coarsen, if needed
            if level == 0:
                zyx_new = data_czyx[ind_chl]
            else:
                zyx_new = da.coarsen(
                    np.mean,
                    pyramid[level - 1][ind_chl],
                    {1: coarsening_xy, 2: coarsening_xy},
                    trim_excess=True,
                )
            list_data_channels.append(zyx_new)

        # Stack several channels
        data_czyx_new = da.stack(list_data_channels, axis=0)

        # Rechunk if needed
        if apply_rechunking:
            data_czyx_final = data_czyx_new.rechunk(
                {2: chunk_size_y, 3: chunk_size_x}, balance=True
            )
        else:
            data_czyx_final = data_czyx_new

        # Append new level
        pyramid.append(data_czyx_final)

    return pyramid
