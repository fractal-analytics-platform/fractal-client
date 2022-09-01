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
from pathlib import Path
from typing import Any
from typing import Dict
from typing import Iterable
from typing import Optional

import dask.array as da
from devtools import debug

from fractal.tasks.lib_pyramid_creation import write_pyramid


def maximum_intensity_projection(
    *,
    input_paths: Iterable[Path],
    output_path: Path,
    metadata: Optional[Dict[str, Any]] = None,
    component: str = None,
):

    """
    Perform maximum-intensity projection along Z axis, and store the output in
    a new zarr file.

    Examples
      input_paths[0] = /tmp/out_mip/*.zarr  (Path)
      output_path = /tmp/out_mip/*zarr   (Path)
      metadata = {"num_levels": 2, "coarsening_xy": 2, ...}
      component = plate.zarr/B/03/0     (str)
    """

    # Preliminary checks
    if len(input_paths) > 1:
        raise NotImplementedError

    # Read some parameters from metadata
    num_levels = metadata["num_levels"]
    coarsening_xy = metadata["coarsening_xy"]
    plate, well = component.split(".zarr/")

    zarrurl_old = metadata["replicate_zarr"]["sources"][plate] + "/" + well
    clean_output_path = output_path.parent.resolve()
    zarrurl_new = (clean_output_path / component).as_posix()
    debug(zarrurl_old)
    debug(zarrurl_new)

    # Hard-coded values (by now) of chunk sizes to be passed to rechunk,
    # both at level 0 (before coarsening) and at levels 1,2,.. (after
    # repeated coarsening).
    # Note that balance=True may override these values.
    # FIXME should this be inferred from somewhere?
    chunk_size_x = 2560
    chunk_size_y = 2160

    # Load 0-th level
    data_czyx = da.from_zarr(zarrurl_old + "/0")
    num_channels = data_czyx.shape[0]
    # Loop over channels
    accumulate_chl = []
    for ind_ch in range(num_channels):
        # Perform MIP for each channel of level 0
        mip_yx = da.stack([da.max(data_czyx[ind_ch], axis=0)], axis=0)
        accumulate_chl.append(mip_yx)
    accumulate_chl = da.stack(accumulate_chl, axis=0)

    # Construct resolution pyramid
    write_pyramid(
        accumulate_chl,
        newzarrurl=zarrurl_new,
        overwrite=False,
        coarsening_xy=coarsening_xy,
        num_levels=num_levels,
        chunk_size_x=chunk_size_x,
        chunk_size_y=chunk_size_y,
    )


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(prog="maximum_intensity_projection.py")
    parser.add_argument(
        "-z", "--zarrurl", help="zarr url, at the FOV level", required=True
    )

    parser.add_argument(
        "-cxy",
        "--coarsening_xy",
        default=2,
        type=int,
        help="coarsening factor along X and Y (optional, defaults to 2)",
    )

    args = parser.parse_args()
    maximum_intensity_projection(args.zarrurl, args.coarsening_xy)
