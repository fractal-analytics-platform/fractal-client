import math
from typing import Dict
from typing import List
from typing import Union

import anndata as ad
import numpy as np
import pandas as pd


def prepare_ROIs_table(
    df: pd.DataFrame, image_size: Union[Dict, None] = None
) -> ad.AnnData:
    if image_size is None:
        raise Exception("Missing image_size arg in prepare_ROIs_table")

    df["x_micrometer"] -= df["x_micrometer"].min()
    df["y_micrometer"] -= df["y_micrometer"].min()
    df["z_micrometer"] -= df["z_micrometer"].min()

    df["len_x_micrometer"] = image_size["x"] * df["pixel_size_x"]
    df["len_y_micrometer"] = image_size["y"] * df["pixel_size_y"]
    df["len_z_micrometer"] = df["pixel_size_z"]

    df.drop("bit_depth", inplace=True, axis=1)

    df = df.astype(np.float32)

    adata = ad.AnnData(X=df, dtype=np.float32)
    adata.obs_names = [f"FOV_{i+1:d}" for i in range(adata.n_obs)]
    adata.var_names = df.columns

    return adata


def convert_ROI_table_to_indices(
    ROI: ad.AnnData, level: int = 0, coarsening_xy: int = 2
) -> List[List]:
    list_indices = []

    for FOV in range(ROI.n_obs):

        # Extract data from anndata table
        # FIXME: is there a better way to do this??
        x_micrometer = ROI[FOV, ROI.var_names == "x_micrometer"].X[0, 0]
        y_micrometer = ROI[FOV, ROI.var_names == "y_micrometer"].X[0, 0]
        z_micrometer = ROI[FOV, ROI.var_names == "z_micrometer"].X[0, 0]
        len_x_micrometer = ROI[FOV, ROI.var_names == "len_x_micrometer"].X[
            0, 0
        ]
        len_y_micrometer = ROI[FOV, ROI.var_names == "len_y_micrometer"].X[
            0, 0
        ]
        len_z_micrometer = ROI[FOV, ROI.var_names == "len_z_micrometer"].X[
            0, 0
        ]
        pixel_size_x = ROI[FOV, ROI.var_names == "pixel_size_x"].X[0, 0]
        pixel_size_y = ROI[FOV, ROI.var_names == "pixel_size_y"].X[0, 0]
        pixel_size_z = ROI[FOV, ROI.var_names == "pixel_size_z"].X[0, 0]

        # Set pyramid-level pixel sizes
        prefactor = coarsening_xy**level
        pixel_size_x *= prefactor
        pixel_size_y *= prefactor

        # Identify indices along the three dimensions
        # FIXME: We don't need Z indices for FOVs, but perhaps we will for
        #        more complex 3D shapes
        start_x = x_micrometer / pixel_size_x
        end_x = (x_micrometer + len_x_micrometer) / pixel_size_x
        start_y = y_micrometer / pixel_size_y
        end_y = (y_micrometer + len_y_micrometer) / pixel_size_y
        start_z = z_micrometer / pixel_size_z
        end_z = (z_micrometer + len_z_micrometer) / pixel_size_z
        indices = [start_z, end_z, start_y, end_y, start_x, end_x]

        # Round indices to lower integer
        # FIXME: to be checked
        indices = list(map(math.floor, indices))

        list_indices.append(indices)

    return list_indices
