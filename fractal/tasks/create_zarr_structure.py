import os
from glob import glob

import zarr

from fractal.tasks.lib_parse_filename_metadata import parse_metadata


def create_zarr_structure(
    in_path=None,
    out_path=None,
    ext=None,
    num_levels=None,
):

    """
    Create (and store) the zarr folder, without reading or writing data.


    :param in_path: path of images
    :type in_path: str
    :param out_path: path for output zarr files
    :type out_path: str
    :param ext: extension of images (e.g. tiff, png, ..)
    :type ext: str
    :param num_levels: number of coarsening levels in the pyramid
    :type num_levels: int
    """

    # Find all plates and all channels
    plates = []
    channels = []
    if not in_path.endswith("/"):
        in_path += "/"
    for i in glob(in_path + "*." + ext):
        try:
            metadata = parse_metadata(os.path.basename(i))
            plates.append(metadata["plate"])
            channels.append(f"A{metadata['A']}_C{metadata['C']}")
        except IndexError:
            print("IndexError for ", i)
            pass
    plates = sorted(list(set(plates)))
    channels = sorted(list(set(channels)))
    print("Find all plates/channels in", in_path + "*." + ext)
    print(f"Plates:   {plates}")
    print(f"Channels: {channels}")

    # HARDCODED CHANNELS AND THEIR PROPERTIES
    allowed_channels = ["A01_C01", "A01_C02", "A02_C03"]
    labels_allowed_channel = {
        "A01_C01": "DAPI",
        "A01_C02": "nanog",
        "A02_C03": "Lamin B1",
    }
    colormaps = ["00FFFF", "FF00FF", "FFFF00"]

    # Check that all channels are in the allowed_channels
    if not set(channels).issubset(set(allowed_channels)):
        msg = "ERROR in create_zarr_structure\n"
        msg += f"channels: {channels}\n"
        msg += f"allowed_channels: {allowed_channels}\n"
        raise Exception(msg)

    # Sort channels according to allowed_channels, and assign increasing index
    # actual_channels is a list of entries like (0, "A01_C01", "DAPI")
    actual_channels = []
    ind_channel = 0
    print("-" * 80)
    print("actual_channels")
    print("-" * 80)
    for ch in allowed_channels:
        if ch in channels:
            actual_channels.append(
                [ind_channel, ch, labels_allowed_channel[ch]]
            )
            ind_channel += 1
            print(ind_channel, ch, labels_allowed_channel[ch])
    print("-" * 80)

    well = []

    zarrurls = {"plate": [], "well": []}

    # loop over plate, each plate could have n wells
    # debug(plate_unique)
    for plate in plates:
        group_plate = zarr.group(out_path + f"{plate}.zarr")
        zarrurls["plate"].append(out_path + f"{plate}.zarr")
        well = [
            parse_metadata(os.path.basename(fn))["well"]
            for fn in glob(in_path + f"{plate}_*." + ext)
        ]
        well_unique = set(well)

        well_rows_columns = [
            ind for ind in sorted([(n[0], n[1:]) for n in well_unique])
        ]

        group_plate.attrs["plate"] = {
            "acquisitions": [
                {"id": id_, "name": name} for id_, name in enumerate(plates)
            ],
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

            # Assumption: All channels in actual_channels are present #FIXME

            group_well.attrs["well"] = {
                "images": [
                    {"path": "0"}  # multiscale level, until pyramids just 0
                ],
                "version": "0.3",
            }

            group_field = group_well.create_group("0/")  # noqa: F841
            zarrurls["well"].append(
                out_path + f"{plate}.zarr/{row}/{column}/0/"
            )

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
                    "datasets": [
                        {"path": f"{level}"} for level in range(num_levels)
                    ],
                }
            ]
            group_field.attrs["omero"] = {
                "id": 1,
                "name": "TBD",
                "version": "0.4",
                "channels": [
                    {
                        # "active": "rue,
                        "coefficient": 1,
                        "color": colormaps[channel[0]],
                        "family": "linear",
                        # "inverted": false,
                        "label": channel[2],
                        "window": {
                            "end": 600,
                            "max": 65535,
                            "min": 0,
                            "start": 110,
                        },
                    }
                    for channel in actual_channels
                ],
            }

    return zarrurls, actual_channels


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(prog="create_zarr_structure")
    parser.add_argument(
        "-i", "--in_path", help="directory containing the input files"
    )
    parser.add_argument(
        "-o", "--out_path", help="directory for the outnput zarr files"
    )
    parser.add_argument(
        "-e",
        "--ext",
        help="source images extension",
    )
    parser.add_argument(
        "-nl",
        "--num_levels",
        type=int,
        help="number of levels in the Zarr pyramid",
    )

    args = parser.parse_args()
    create_zarr_structure(
        in_path=args.in_path,
        out_path=args.out_path,
        ext=args.ext,
        num_levels=args.num_levels,
    )
