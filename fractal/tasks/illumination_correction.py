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
import json
import warnings
from concurrent.futures import ThreadPoolExecutor

import anndata as ad
import dask
import dask.array as da
import numpy as np
from skimage.io import imread

from fractal.tasks.lib_pyramid_creation import write_pyramid
from fractal.tasks.lib_regions_of_interest import convert_ROI_table_to_indices
from fractal.tasks.lib_regions_of_interest import (
    split_3D_indices_into_z_layers,
)
from fractal.tasks.lib_zattrs_utils import extract_zyx_pixel_sizes


def correct(
    img,
    illum_img=None,
    background=110,
):
    """
    Corrects single Z level input image using an illumination profile

    The illumination profile image can be uint8 or uint16.
    It needs to follow the illumination profile. e.g. bright in the
    center of the image, dim outside

    :param img: input image to be corrected, either uint8 or uint16
    :type img: np.array
    :param illum_img: correction matrix
    :type illum_img: np.array
    :param background: value for background subtraction (optional, default 110)
    :type background: int

    """

    # Check shapes
    if illum_img.shape != img.shape[1:]:
        raise Exception(
            "Error in illumination_correction\n"
            f"img.shape: {img.shape}\n"
            f"illum_img.shape: {illum_img.shape}"
        )

    # Background subtraction
    # FIXME: is there a problem with these changes?
    # devdoc.net/python/dask-2.23.0-doc/delayed-best-practices.html
    # ?highlight=delayed#don-t-mutate-inputs
    img[img <= background] = 0
    img[img > background] = img[img > background] - background

    # Apply the illumination correction
    # (normalized by the max value in the illum_img)
    img_corr = img / (illum_img / np.max(illum_img))[None, :, :]

    # Handle edge case: The illumination correction can increase a value
    # beyond the limit of the encoding, e.g. beyond 65535 for 16bit
    # images. This clips values that surpass this limit and triggers
    # a warning
    if np.sum(img_corr > np.iinfo(img.dtype).max) > 0:
        warnings.warn(
            f"The illumination correction created values \
                       beyond the max range of your current image \
                       type. These have been clipped to \
                       {np.iinfo(img.dtype).max}"
        )
        img_corr[img_corr > np.iinfo(img.dtype).max] = np.iinfo(img.dtype).max

    return img_corr.astype(img.dtype)


def illumination_correction(
    zarrurl,
    overwrite=False,
    newzarrurl=None,
    chl_list=None,
    path_dict_corr=None,
    coarsening_xy=2,
    background=110,
    num_threads=10,
):

    """
    Perform illumination correction of the array in zarrurl.

    :param zarrurl: input zarr, at the site level (e.g. x.zarr/B/03/0/)
    :type zarrurl: str
    :param overwrite: whether to overwrite an existing zarr file
    :type overwrite: bool
    :param newzarrurl: output zarr, at the site level (e.g. x.zarr/B/03/0/)
    :type newzarrurl: str
    :param chl_list: list of channel names (e.g. A01_C01)
    :type chl_list: list
    :param path_dict_corr: path of JSON file with info on illumination matrices
    :type path_dict_corr: str
    :param coarsening_xy: coarsening factor in XY (optional, default 2)
    :type coarsening_xy: xy
    :param background: value for background subtraction (optional, default 110)
    :type background: int

    """

    # Check that only one output option is chosen
    if overwrite and (newzarrurl is not None) and (newzarrurl != zarrurl):
        raise Exception(
            "ERROR in illumination_correction: "
            f"overwrite={overwrite} and newzarrurl={newzarrurl}."
        )
    if newzarrurl is None and not overwrite:
        raise Exception(
            "ERROR in illumination_correction: "
            f"overwrite={overwrite} and newzarrurl={newzarrurl}."
        )

    # Sanitize zarr paths
    if not zarrurl.endswith("/"):
        zarrurl += "/"
    if overwrite:
        newzarrurl = zarrurl
    else:
        if not newzarrurl.endswith("/"):
            newzarrurl += "/"

    # Read FOV ROIs
    FOV_ROI_table = ad.read_zarr(f"{zarrurl}tables/FOV_ROI_table")

    # Read pixel sizes from zattrs file
    full_res_pxl_sizes_zyx = extract_zyx_pixel_sizes(
        zarrurl + ".zattrs", level=0
    )

    # Create list of indices for 3D FOVs spanning the entire Z direction
    list_indices = convert_ROI_table_to_indices(
        FOV_ROI_table,
        level=0,
        coarsening_xy=coarsening_xy,
        full_res_pxl_sizes_zyx=full_res_pxl_sizes_zyx,
    )

    # Extract image size from FOV-ROI indices
    # Note: this works at level=0, where FOVs should all be of the exact same
    #       size (in pixels)
    ref_img_size = None
    for indices in list_indices:
        img_size = (indices[3] - indices[2], indices[5] - indices[4])
        if ref_img_size is None:
            ref_img_size = img_size
        else:
            if img_size != ref_img_size:
                raise Exception(
                    "ERROR: inconsistent image sizes in list_indices"
                )
    img_size_y, img_size_x = img_size[:]

    # Load paths of correction matrices
    with open(path_dict_corr, "r") as jsonfile:
        dict_corr = json.load(jsonfile)
    root_path_corr = dict_corr.pop("root_path_corr")
    if not root_path_corr.endswith("/"):
        root_path_corr += "/"

    # Assemble dictionary of matrices and check their shapes
    corrections = {}
    for ind_ch, ch in enumerate(chl_list):
        corrections[ch] = imread(root_path_corr + dict_corr[ch])
        if corrections[ch].shape != (img_size_y, img_size_x):
            raise Exception(
                "Error in illumination_correction, "
                "correction matrix has wrong shape."
            )

    # Read number of levels from .zattrs of original zarr file
    with open(zarrurl + ".zattrs", "r") as inputjson:
        zattrs = json.load(inputjson)
    num_levels = len(zattrs["multiscales"][0]["datasets"])

    # Load highest-resolution level from original zarr array
    data_czyx = da.from_zarr(zarrurl + "/0")
    dtype = data_czyx.dtype

    # Check that input array is made of images (in terms of shape/chunks)
    nc, nz, ny, nx = data_czyx.shape
    if (ny % img_size_y != 0) or (nx % img_size_x != 0):
        raise Exception(
            "Error in illumination_correction, "
            f"data_czyx.shape: {data_czyx.shape}"
        )
    chunks_c, chunks_z, chunks_y, chunks_x = data_czyx.chunks
    if len(set(chunks_c)) != 1 or chunks_c[0] != 1:
        raise Exception(
            f"Error in illumination_correction, chunks_c: {chunks_c}"
        )
    if len(set(chunks_z)) != 1 or chunks_z[0] != 1:
        raise Exception(
            f"Error in illumination_correction, chunks_z: {chunks_z}"
        )
    if len(set(chunks_y)) != 1 or chunks_y[0] != img_size_y:
        raise Exception(
            f"Error in illumination_correction, chunks_y: {chunks_y}"
        )
    if len(set(chunks_x)) != 1 or chunks_x[0] != img_size_x:
        raise Exception(
            f"Error in illumination_correction, chunks_x: {chunks_x}"
        )

    # Create the final list of single-Z-layer FOVs
    list_indices = split_3D_indices_into_z_layers(list_indices)

    # Prepare delayed function
    delayed_correct = dask.delayed(correct)

    # Loop over channels
    data_czyx_new = []
    for ind_ch, ch in enumerate(chl_list):
        # Set correction matrix
        illum_img = corrections[ch]
        # Select data for given channel
        data_zyx = data_czyx[ind_ch]
        # Initialize empty array for corrected images
        data_zyx_new = da.empty_like(data_zyx)
        # Loop over FOVs
        for indices in list_indices:
            s_z, e_z, s_y, e_y, s_x, e_x = indices[:]
            shape = [e_z - s_z, e_y - s_y, e_x - s_x]
            if min(shape) == 0:
                raise Exception(f"ERROR: ROI indices lead to shape {shape}")
            new_img = delayed_correct(
                data_zyx[s_z:e_z, s_y:e_y, s_x:e_x],
                illum_img,
                background=background,
            )
            data_zyx_new[s_z:e_z, s_y:e_y, s_x:e_x] = da.from_delayed(
                new_img, shape, dtype
            )
        data_czyx_new.append(data_zyx_new)

    # Combine all 3D channel arrays into a single 4D array
    accumulated_data = da.stack(data_czyx_new, axis=0)

    # Construct resolution pyramid
    with dask.config.set(pool=ThreadPoolExecutor(num_threads)):
        write_pyramid(
            accumulated_data,
            newzarrurl=newzarrurl,
            overwrite=overwrite,
            coarsening_xy=coarsening_xy,
            num_levels=num_levels,
            chunk_size_x=img_size_x,
            chunk_size_y=img_size_y,
        )


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(prog="illumination_correction.py")
    parser.add_argument(
        "-z", "--zarrurl", help="zarr url, at the FOV level", required=True
    )
    parser.add_argument(
        "-o",
        "--overwrite",
        action="store_true",
        help="overwrite original zarr file",
    )
    parser.add_argument(
        "-znew",
        "--newzarrurl",
        help="path of the new zarr file",
    )
    parser.add_argument(
        "--path_dict_corr",
        help="path of JSON file with info on illumination matrices",
    )
    parser.add_argument(
        "-C",
        "--chl_list",
        nargs="+",
        help="list of channel names (e.g. A01_C01)",
    )
    parser.add_argument(
        "-cxy",
        "--coarsening_xy",
        default=2,
        type=int,
        help="coarsening factor along X and Y (optional, defaults to 2)",
    )
    parser.add_argument(
        "-bg",
        "--background",
        default=110,
        type=int,
        help=(
            "threshold for background subtraction"
            " (optional, defaults to 110)"
        ),
    )

    args = parser.parse_args()
    illumination_correction(
        args.zarrurl,
        overwrite=args.overwrite,
        newzarrurl=args.newzarrurl,
        path_dict_corr=args.path_dict_corr,
        chl_list=args.chl_list,
        coarsening_xy=args.coarsening_xy,
        background=args.background,
    )
