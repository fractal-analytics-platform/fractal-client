import math
from typing import Dict
from typing import List
from typing import Union

import anndata as ad
import numpy as np
import pandas as pd


def prepare_FOV_ROI_table(
    df: pd.DataFrame,
    well_size: Union[Dict, None] = None,
) -> ad.AnnData:
    if well_size is None:
        raise Exception("Missing well_size arg in prepare_ROIs_table")

    # Reset reference values for coordinates
    df["x_micrometer"] -= df["x_micrometer"].min()
    df["y_micrometer"] -= df["y_micrometer"].min()
    df["z_micrometer"] -= df["z_micrometer"].min()

    # Obtain box size in physical units
    df["len_x_micrometer"] = well_size["x"] * df["pixel_size_x"]
    df["len_y_micrometer"] = well_size["y"] * df["pixel_size_y"]
    df["len_z_micrometer"] = well_size["z"] * df["pixel_size_z"]

    # Remove unused column
    df.drop("bit_depth", inplace=True, axis=1)

    # Assign dtype explicitly, to avoid
    # >> UserWarning: X converted to numpy array with dtype float64
    # when creating AnnData object
    df = df.astype(np.float32)

    # Convert DataFrame index to str, to avoid
    # >> ImplicitModificationWarning: Transforming to str index
    # when creating AnnData object
    df.index = df.index.astype(str)

    # Create an AnnData object directly from the DataFrame
    adata = ad.AnnData(X=df)

    # Rename rows and columns
    adata.obs_names = [f"FOV_{i+1:d}" for i in range(adata.n_obs)]
    adata.var_names = list(map(str, df.columns))

    return adata


def convert_FOV_ROIs_3D_to_2D(adata: ad.AnnData = None) -> ad.AnnData:

    # FIXME: use this piece of code (or similar) in tests
    """
    adata = ad.AnnData(X=np.zeros((2, 3)), dtype=int)
    adata.obs_names = ["FOV1", "FOV2"]
    adata.var_names = ["z", "len_z", "pixel_size_z"]
    adata[:, "z"].X = [[2], [4]]
    adata[:, "len_z"].X = [[5], [7]]
    adata[:, "pixel_size_z"].X = [[2], [2]]
    """

    # Compress a 3D stack of images to a single Z plane,
    # with thickness equal to pixel_size_z
    df = adata.to_df()
    df["len_z"] = df["pixel_size_z"]

    # Assign dtype explicitly, to avoid
    # >> UserWarning: X converted to numpy array with dtype float64
    # when creating AnnData object
    df = df.astype(np.float32)

    # Create an AnnData object directly from the DataFrame
    new_adata = ad.AnnData(X=df)

    # Rename rows and columns
    new_adata.obs_names = [f"FOV_{i+1:d}" for i in range(new_adata.n_obs)]
    new_adata.var_names = list(map(str, df.columns))

    return new_adata


def convert_ROI_table_to_indices(
    ROI: ad.AnnData,
    level: int = 0,
    coarsening_xy: int = 2,
) -> List[List]:

    list_indices = []

    for FOV in sorted(ROI.obs_names):

        # Extract data from anndata table
        x_micrometer = ROI[FOV, "x_micrometer"].X[0, 0]
        y_micrometer = ROI[FOV, "y_micrometer"].X[0, 0]
        z_micrometer = ROI[FOV, "z_micrometer"].X[0, 0]
        len_x_micrometer = ROI[FOV, "len_x_micrometer"].X[0, 0]
        len_y_micrometer = ROI[FOV, "len_y_micrometer"].X[0, 0]
        len_z_micrometer = ROI[FOV, "len_z_micrometer"].X[0, 0]
        pixel_size_x = ROI[FOV, "pixel_size_x"].X[0, 0]
        pixel_size_y = ROI[FOV, "pixel_size_y"].X[0, 0]
        pixel_size_z = ROI[FOV, "pixel_size_z"].X[0, 0]

        # Set pyramid-level pixel sizes
        prefactor = coarsening_xy**level
        pixel_size_x *= prefactor
        pixel_size_y *= prefactor

        # Identify indices along the three dimensions
        start_x = x_micrometer / pixel_size_x
        end_x = (x_micrometer + len_x_micrometer) / pixel_size_x
        start_y = y_micrometer / pixel_size_y
        end_y = (y_micrometer + len_y_micrometer) / pixel_size_y
        start_z = z_micrometer / pixel_size_z
        end_z = (z_micrometer + len_z_micrometer) / pixel_size_z
        indices = [start_z, end_z, start_y, end_y, start_x, end_x]

        # Round indices to lower integer
        # FIXME: to be checked/tested
        indices = list(map(math.floor, indices))

        # Append ROI indices to to list
        list_indices.append(indices)

    return list_indices
