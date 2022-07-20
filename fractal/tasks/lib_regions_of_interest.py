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
    ROI: ad.AnnData,
    level: int = 0,
    coarsening_xy: int = 2,
    num_z_replicas: int = 1,
) -> List[List]:
    list_indices = []

    for FOV in ROI.obs_names:

        # Extract data from anndata table
        # FIXME: is there a better way to do this??
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

        # Default behavior
        if num_z_replicas == 1:
            list_indices.append(indices)
        # Create 3D stack of 2D ROIs
        else:
            # Check that this ROI is 2D, i.e. it has z indices [0:1]
            if start_z != 0 or end_z != 1:
                raise Exception(
                    f"ERROR: num_z_replicas={num_z_replicas}, "
                    f"but [start_z,end_z]={[start_z,end_z]}"
                )
            # Loop over Z planes
            for z_start in range(num_z_replicas):
                indices[0:2] = [z_start, z_start + 1]
                list_indices.append(indices[:])

    return list_indices
