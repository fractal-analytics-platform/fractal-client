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
from glob import glob

import zarr


def replicate_zarr_structure(zarrurl, newzarrurl=None):
    """
    Duplicate an input zarr structure to a new path.

    :param zarrurl: structure of the input zarr folder
    :type zarrurl: str
    :param newzarrurl: structure of the input zarr folder
    :type newzarrurl: str
    """

    # Sanitize and check input zarr path
    if not zarrurl.endswith("/"):
        zarrurl += "/"
    if not zarrurl.endswith(".zarr/"):
        raise Exception(
            "Error in replicate_zarr_structure, "
            f"zarrurl={zarrurl} does not end with .zarr/"
        )

    # Sanitize and check output zarr path
    if newzarrurl is None:
        newzarrurl = zarrurl[-5:] + "_corr.zarr/"
    else:
        if not newzarrurl.endswith("/"):
            newzarrurl += "/"
        if not newzarrurl.endswith(".zarr/"):
            raise Exception(
                "Error in replicate_zarr_structure, "
                f"newzarrurl={newzarrurl} does not end with .zarr/"
            )

    # Identify properties of input zarr file
    well_rows_columns = sorted(
        [rc.split("/")[-2:] for rc in glob(zarrurl + "*/*")]
    )
    levels = sorted(
        list(set([rc.split("/")[-1] for rc in glob(zarrurl + "*/*/*/*")]))
    )

    group_plate = zarr.group(newzarrurl)
    plate = zarrurl.replace(".zarr/", "").split("/")[-1]
    group_plate.attrs["plate"] = {
        "acquisitions": [{"id": 0, "name": plate}],
        # takes unique cols from (row,col) tuples
        "columns": sorted(
            [
                {"name": u_col}
                for u_col in set(
                    [
                        well_row_column[1]
                        for well_row_column in well_rows_columns
                    ]
                )
            ],
            key=lambda key: key["name"],
        ),
        # takes unique rows from (row,col) tuples
        "rows": sorted(
            [
                {"name": u_row}
                for u_row in set(
                    [
                        well_row_column[0]
                        for well_row_column in well_rows_columns
                    ]
                )
            ],
            key=lambda key: key["name"],
        ),
        "wells": [
            {
                "path": well_row_column[0] + "/" + well_row_column[1],
            }
            for well_row_column in well_rows_columns
        ],
    }

    for row, column in well_rows_columns:

        group_well = group_plate.create_group(f"{row}/{column}/")

        group_well.attrs["well"] = {
            "images": [{"path": "0"}],  # FOV (only 0, by now)
            "version": "0.3",
        }
        group_field = group_well.create_group("0/")  # noqa: F841

        group_field.attrs["multiscales"] = [
            {
                "version": "0.3",
                "axes": [
                    {"name": "c", "type": "channel"},
                    {
                        "name": "z",
                        "type": "space",
                        "unit": "micrometer",
                    },
                    {"name": "y", "type": "space"},
                    {"name": "x", "type": "space"},
                ],
                "datasets": [{"path": level} for level in levels],
            }
        ]

        # Copy .zattrs file at the COL/ROW/SITE level
        path_zattrs = zarrurl + f"{row}/{column}/0/.zattrs"
        with open(path_zattrs) as zattrs_file:
            zattrs = json.load(zattrs_file)
            group_field.attrs["omero"] = zattrs["omero"]


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(prog="FIXME")
    parser.add_argument(
        "-zo", "--zarrurl_old", help="old zarr url", required=True
    )
    parser.add_argument(
        "-zn", "--zarrurl_new", help="new zarr url", required=True
    )
    args = parser.parse_args()
    replicate_zarr_structure(args.zarrurl_old, newzarrurl=args.zarrurl_new)
