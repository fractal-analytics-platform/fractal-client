from glob import glob

import zarr


def replicate_zarr_structure_mip(
    zarrurl,
):

    """
    FIXME
    input: file.zarr
    """

    if not zarrurl.endswith("/"):
        zarrurl += "/"

    if not zarrurl.endswith(".zarr/"):
        raise

    zarrurl_mip = zarrurl.replace(".zarr/", "_mip.zarr/")
    plate = zarrurl.replace(".zarr/", "").split("/")[-1]

    group_plate = zarr.group(zarrurl_mip)

    well_rows_columns = sorted(
        [rc.split("/")[-2:] for rc in glob(zarrurl + "*/*")]
    )
    levels = sorted(
        list(set([rc.split("/")[-1] for rc in glob(zarrurl + "*/*/*/*")]))
    )

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
    # debug(well_rows_columns)

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


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(prog="FIXME")
    parser.add_argument("-z", "--zarrurl", help="zarr url", required=True)
    args = parser.parse_args()
    replicate_zarr_structure_mip(args.zarrurl)
