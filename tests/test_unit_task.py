import json

from devtools import debug

from fractal.tasks import __FRACTAL_MANIFEST__
from fractal.tasks.create_zarr_structure import create_zarr_structure

create_zarr_structure_manifest = next(
    item
    for item in __FRACTAL_MANIFEST__
    if item["name"] == "Create OME-ZARR structure"
)


def test_create_zarr_structure(tmp_path, testdata_path):
    input_paths = [testdata_path / "png/*.png"]
    output_path = tmp_path
    default_args = create_zarr_structure_manifest["default_args"]

    from glob import glob

    debug(glob(input_paths[0].as_posix()))

    zarrurls, actual_channels = create_zarr_structure(
        input_paths=input_paths, output_path=output_path, **default_args
    )

    debug(list(output_path.glob("*")))
    zattrs = output_path / "myplate.zarr/.zattrs"
    with open(zattrs) as f:
        data = json.load(f)
        debug(data)
    assert len(data["plate"]["wells"]) == 1

    debug(zarrurls["well"])
