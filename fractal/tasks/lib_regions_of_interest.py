import json
import math
from typing import List
from typing import Tuple
from typing import Union

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
    pixel_sizes_zyx: Union[List[float], Tuple[float]] = None,
) -> List[List[int]]:

    list_indices = []

    pixel_size_z, pixel_size_y, pixel_size_x = pixel_sizes_zyx

    for FOV in sorted(ROI.obs_names):

        # Extract data from anndata table
        x_micrometer = ROI[FOV, "x_micrometer"].X[0, 0]
        y_micrometer = ROI[FOV, "y_micrometer"].X[0, 0]
        z_micrometer = ROI[FOV, "z_micrometer"].X[0, 0]
        len_x_micrometer = ROI[FOV, "len_x_micrometer"].X[0, 0]
        len_y_micrometer = ROI[FOV, "len_y_micrometer"].X[0, 0]
        len_z_micrometer = ROI[FOV, "len_z_micrometer"].X[0, 0]

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


def split_3D_indices_into_z_layers(
    list_indices: List[List[int]],
) -> List[List[int]]:

    num_z_layers = None
    new_list_indices = []
    for indices in list_indices:
        if num_z_layers is None:
            num_z_layers = indices[1]
        else:
            if indices[1] != num_z_layers:
                raise Exception(
                    "Inconsistent num_z_layers in split_indices_into_2D_layers"
                )
        for ind_z in range(num_z_layers):
            new_indices = [ind_z, ind_z + 1] + indices[2:]
            new_list_indices.append(new_indices)

    return new_list_indices


def _inspect_ROI_table(
    path: str = None,
    level: int = 0,
    coarsening_xy: int = 2,
    pixel_sizes_zyx=[1.0, 0.1625, 0.1625],
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
        pixel_sizes_zyx=pixel_sizes_zyx,
    )

    list_indices = split_3D_indices_into_z_layers(list_indices)

    print(f"level:         {level}")
    print(f"coarsening_xy: {coarsening_xy}")
    print("list_indices:")
    for indices in list_indices:
        print(indices)
    print()


def temporary_test():

    pixel_size_z = 1.0
    pixel_size_y = 0.1625
    pixel_size_x = 0.1625

    df = pd.DataFrame(np.zeros((2, 10)), dtype=int)
    df.index = ["FOV1", "FOV2"]
    df.columns = [
        "x_micrometer",
        "y_micrometer",
        "z_micrometer",
        "x_pixel",
        "y_pixel",
        "z_pixel",
        "pixel_size_x",
        "pixel_size_y",
        "pixel_size_z",
        "bit_depth",
    ]
    df["x_micrometer"] = [0.0, 416.0]
    df["y_micrometer"] = [0.0, 0.0]
    df["z_micrometer"] = [0.0, 0.0]
    df["x_pixel"] = [2560, 2560]
    df["y_pixel"] = [2160, 2160]
    df["z_pixel"] = [5, 5]
    df["pixel_size_x"] = [pixel_size_x] * 2
    df["pixel_size_y"] = [pixel_size_y] * 2
    df["pixel_size_z"] = [pixel_size_z] * 2
    df["bit_depth"] = [16.0, 16.0]

    print("DataFrame")
    print(df)
    print()

    adata = prepare_FOV_ROI_table(df)

    print("AnnData table")
    print(adata.var_names)
    print(adata.obs_names)
    print(adata.X)
    print()

    print("Indices 3D")
    pixel_sizes_zyx = [pixel_size_z, pixel_size_y, pixel_size_x]
    list_indices = convert_ROI_table_to_indices(
        adata, level=0, coarsening_xy=2, pixel_sizes_zyx=pixel_sizes_zyx
    )
    for indices in list_indices:
        print(indices)
    print()

    print("Indices 3D / split")
    list_indices = split_3D_indices_into_z_layers(list_indices)
    for indices in list_indices:
        print(indices)
    print()

    print("Indices 2D")
    adata = convert_FOV_ROIs_3D_to_2D(adata, pixel_size_z)
    list_indices = convert_ROI_table_to_indices(
        adata, level=0, coarsening_xy=2, pixel_sizes_zyx=pixel_sizes_zyx
    )
    for indices in list_indices:
        print(indices)
    print()


def extract_zyx_pixel_sizes_from_zattrs(zattrs_path: str, level: int = 0):
    with open(zattrs_path, "r") as jsonfile:
        zattrs = json.load(jsonfile)

    try:

        # Identify multiscales
        multiscales = zattrs["multiscales"]

        # Check that there is a single multiscale
        if len(multiscales) > 1:
            raise Exception(f"ERROR: There are {len(multiscales)} multiscales")

        # Check that there are no datasets-global transformations
        if "coordinateTransformations" in multiscales[0].keys():
            raise Exception(
                "ERROR: coordinateTransformations at the multiscales "
                "level are not currently supported"
            )

        # Identify all datasets (AKA pyramid levels)
        datasets = multiscales[0]["datasets"]

        # Select highest-resolution dataset
        transformations = datasets[level]["coordinateTransformations"]
        for t in transformations:
            if t["type"] == "scale":
                return t["scale"]
        raise Exception(
            "ERROR:"
            f" no scale transformation found for level {level}"
            f" in {zattrs_path}"
        )

    except KeyError as e:
        raise KeyError(
            "extract_zyx_pixel_sizes_from_zattrs failed, for {zattrs_path}\n",
            e,
        )


if __name__ == "__main__":
    # import sys
    # args = sys.argv[1:]
    # types = [str, int, int]
    # args = [types[ix](x) for ix, x in enumerate(args)]
    # _inspect_ROI_table(*args)

    temporary_test()
