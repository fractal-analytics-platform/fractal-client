import numpy as np
import pandas as pd
import pytest

from fractal.tasks.lib_regions_of_interest import convert_FOV_ROIs_3D_to_2D
from fractal.tasks.lib_regions_of_interest import convert_ROI_table_to_indices
from fractal.tasks.lib_regions_of_interest import prepare_FOV_ROI_table


PIXEL_SIZE_X = 0.1625
PIXEL_SIZE_Y = 0.1625
PIXEL_SIZE_Z = 1.0

IMG_SIZE_X = 2560
IMG_SIZE_Y = 2160
NUM_Z_PLANES = 4

FOV_IDS = ["1", "2", "7", "9"]
FOV_NAMES = [f"FOV_{ID}" for ID in FOV_IDS]


def get_metadata_dataframe():
    """
    Create artificial metadata dataframe
    """
    df = pd.DataFrame(np.zeros((4, 11)), dtype=int)
    df.index = FOV_IDS
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
        "time",
    ]
    img_size_x_micrometer = IMG_SIZE_X * PIXEL_SIZE_X
    img_size_y_micrometer = IMG_SIZE_Y * PIXEL_SIZE_Y
    df["x_micrometer"] = [
        0.0,
        img_size_x_micrometer,
        0.0,
        img_size_x_micrometer,
    ]
    df["y_micrometer"] = [
        0.0,
        0.0,
        img_size_y_micrometer,
        img_size_y_micrometer,
    ]
    df["z_micrometer"] = [0.0, 0.0, 0.0, 0.0]
    df["x_pixel"] = [IMG_SIZE_X] * 4
    df["y_pixel"] = [IMG_SIZE_Y] * 4
    df["z_pixel"] = [NUM_Z_PLANES] * 4
    df["pixel_size_x"] = [PIXEL_SIZE_X] * 4
    df["pixel_size_y"] = [PIXEL_SIZE_Y] * 4
    df["pixel_size_z"] = [PIXEL_SIZE_Z] * 4
    df["bit_depth"] = [16.0] * 4
    df["time"] = "2020-08-12 15:36:36.234000+0000"

    return df


list_pxl_sizes = []
list_pxl_sizes.append([PIXEL_SIZE_Z, PIXEL_SIZE_Y, PIXEL_SIZE_X])
list_pxl_sizes.append([val + 1e-8 for val in list_pxl_sizes[0]])
list_pxl_sizes.append([val - 1e-8 for val in list_pxl_sizes[0]])


list_level_coarsening = [
    (0, 2),
    (1, 2),
    (2, 2),
    (3, 2),
    (0, 3),
    (1, 3),
    (2, 3),
    (3, 3),
    (0, 7),
    (1, 7),
    (2, 7),
]

list_params = []
for pxl_sizes in list_pxl_sizes:
    for level_coarsening in list_level_coarsening:
        level, coarsening = level_coarsening[:]
        list_params.append((level, coarsening, pxl_sizes))


@pytest.mark.xfail(reason="May fail because of rounding issues (#159)")
@pytest.mark.parametrize(
    "level,coarsening_xy,full_res_pxl_sizes_zyx", list_params
)
def test_ROI_indices_3D(level, coarsening_xy, full_res_pxl_sizes_zyx):

    metadata_dataframe = get_metadata_dataframe()
    adata = prepare_FOV_ROI_table(metadata_dataframe)
    assert list(adata.obs_names) == FOV_NAMES

    list_indices = convert_ROI_table_to_indices(
        adata,
        level=level,
        coarsening_xy=coarsening_xy,
        full_res_pxl_sizes_zyx=full_res_pxl_sizes_zyx,
    )
    print()
    original_shape = (
        NUM_Z_PLANES,
        2 * IMG_SIZE_Y,
        2 * IMG_SIZE_X // coarsening_xy**level,
    )
    expected_shape = (
        NUM_Z_PLANES,
        2 * IMG_SIZE_Y // coarsening_xy**level,
        2 * IMG_SIZE_X // coarsening_xy**level,
    )
    print(f"Pixel sizes: {full_res_pxl_sizes_zyx}")
    print(f"Original shape: {original_shape}")
    print(f"coarsening_xy={coarsening_xy}, level={level}")
    print(f"Expected shape: {expected_shape}")
    print("FOV-ROI indices:")
    for indices in list_indices:
        print(indices)
    print()

    assert list_indices[0][5] == list_indices[1][4]
    assert list_indices[0][3] == list_indices[2][2]
    assert (
        abs(
            (list_indices[0][5] - list_indices[0][4])
            - (list_indices[1][5] - list_indices[1][4])
        )
        < coarsening_xy
    )
    assert (
        abs(
            (list_indices[0][3] - list_indices[0][2])
            - (list_indices[1][3] - list_indices[1][2])
        )
        < coarsening_xy
    )
    assert abs(list_indices[1][5] - expected_shape[2]) < coarsening_xy
    assert abs(list_indices[2][3] - expected_shape[1]) < coarsening_xy
    for indices in list_indices:
        assert indices[0] == 0
        assert indices[1] == NUM_Z_PLANES


@pytest.mark.xfail(reason="May fail because of rounding issues (#159)")
@pytest.mark.parametrize(
    "level,coarsening_xy,full_res_pxl_sizes_zyx", list_params
)
def test_ROI_indices_2D(level, coarsening_xy, full_res_pxl_sizes_zyx):

    metadata_dataframe = get_metadata_dataframe()
    adata = prepare_FOV_ROI_table(metadata_dataframe)
    adata = convert_FOV_ROIs_3D_to_2D(adata, PIXEL_SIZE_Z)

    list_indices = convert_ROI_table_to_indices(
        adata,
        level=level,
        coarsening_xy=coarsening_xy,
        full_res_pxl_sizes_zyx=full_res_pxl_sizes_zyx,
    )

    for indices in list_indices:
        assert indices[0] == 0
        assert indices[1] == 1
