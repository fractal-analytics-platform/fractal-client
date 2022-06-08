import dask.array as da
import numpy as np


def create_pyramid(
    data_czyx,
    coarsening_z=1,
    coarsening_xy=1,
    num_levels=1,
    chunk_size_x=None,
    chunk_size_y=None,
    chl_list=[],
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
    :param chl_list: list of the channels
    :type chl_list: list
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
    pyramid = {level: [] for level in range(num_levels)}
    for ind_chl, chl in enumerate(chl_list):
        zyx_arrays = {}
        for level in range(num_levels):
            if level == 0:
                zyx_arrays[level] = data_czyx[ind_chl]
            else:
                zyx_arrays[level] = da.coarsen(
                    np.min,
                    zyx_arrays[level - 1],
                    {1: coarsening_xy, 2: coarsening_xy},
                    trim_excess=True,
                )
            # Rechunk array
            if apply_rechunking:
                zyx_arrays[level] = zyx_arrays[level].rechunk(
                    {1: chunk_size_y, 2: chunk_size_x}, balance=True
                )
            pyramid[level].append(zyx_arrays[level])

    # pyramid[level] has dimensions (channel, Z, Y, X)
    for ind_level in range(num_levels):
        pyramid[ind_level] = da.stack(pyramid[ind_level], axis=0)

    return pyramid
