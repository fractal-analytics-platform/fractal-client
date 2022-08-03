import math
from typing import Iterable
from typing import List

import anndata as ad
import numpy as np
import pandas as pd


def prepare_FOV_ROI_table(df: pd.DataFrame) -> ad.AnnData:

    for mu in ["x", "y", "z"]:

        # Reset reference values for coordinates
        df[f"{mu}_micrometer"] -= df[f"{mu}_micrometer"].min()

        # Obtain box size in physical units
        df[f"len_{mu}_micrometer"] = df[f"{mu}_pixel"] * df[f"pixel_size_{mu}"]

        # Remove information about pixel sizes in physical units
        df.drop(f"pixel_size_{mu}", inplace=True, axis=1)

        # Remove information about array size in pixels
        df.drop(f"{mu}_pixel", inplace=True, axis=1)

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


def convert_FOV_ROIs_3D_to_2D(
    adata: ad.AnnData = None, pixel_size_z: float = None
) -> ad.AnnData:

    if pixel_size_z is None:
        raise Exception("Missing pixel_size_z in convert_FOV_ROIs_3D_to_2D")

    # Compress a 3D stack of images to a single Z plane,
    # with thickness equal to pixel_size_z
    df = adata.to_df()
    df["len_z_micrometer"] = pixel_size_z

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
    full_res_pxl_sizes_zyx: Iterable[float] = None,
) -> List[List[int]]:

    # Set pyramid-level pixel sizes
    pxl_size_z, pxl_size_y, pxl_size_x = full_res_pxl_sizes_zyx
    prefactor = coarsening_xy**level
    pxl_size_x *= prefactor
    pxl_size_y *= prefactor

    list_indices = []
    for FOV in sorted(ROI.obs_names):

        # Extract data from anndata table
        x_micrometer = ROI[FOV, "x_micrometer"].X[0, 0]
        y_micrometer = ROI[FOV, "y_micrometer"].X[0, 0]
        z_micrometer = ROI[FOV, "z_micrometer"].X[0, 0]
        len_x_micrometer = ROI[FOV, "len_x_micrometer"].X[0, 0]
        len_y_micrometer = ROI[FOV, "len_y_micrometer"].X[0, 0]
        len_z_micrometer = ROI[FOV, "len_z_micrometer"].X[0, 0]

        # Identify indices along the three dimensions
        start_x = x_micrometer / pxl_size_x
        end_x = (x_micrometer + len_x_micrometer) / pxl_size_x
        start_y = y_micrometer / pxl_size_y
        end_y = (y_micrometer + len_y_micrometer) / pxl_size_y
        start_z = z_micrometer / pxl_size_z
        end_z = (z_micrometer + len_z_micrometer) / pxl_size_z
        indices = [start_z, end_z, start_y, end_y, start_x, end_x]

        # Round indices to lower integer
        indices = list(map(math.floor, indices))

        # Append ROI indices to to list
        list_indices.append(indices[:])

    return list_indices


def _inspect_ROI_table(
    path: str = None,
    level: int = 0,
    coarsening_xy: int = 2,
    full_res_pxl_sizes_zyx=[1.0, 0.1625, 0.1625],
) -> None:

    adata = ad.read_zarr(path)
    df = adata.to_df()
    print("table")
    print(df)
    print()

    list_indices = convert_ROI_table_to_indices(
        adata,
        level=level,
        coarsening_xy=coarsening_xy,
        full_res_pxl_sizes_zyx=full_res_pxl_sizes_zyx,
    )

    print(f"level:         {level}")
    print(f"coarsening_xy: {coarsening_xy}")
    print("list_indices:")
    for indices in list_indices:
        print(indices)
    print()
