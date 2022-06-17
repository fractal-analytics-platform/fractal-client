import json
import os
from glob import glob

import zarr

from fractal.tasks.lib_parse_filename_metadata import parse_metadata


def create_zarr_structure(
    in_path=None,
    out_path=None,
    ext=None,
    path_dict_channels=None,
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

    # FIXME: in_path should be a list!

    # Find all plates and all channels
    plates = []
    channels = []
    dict_plate_prefixes = {}
    if not in_path.endswith("/"):
        in_path += "/"
    for fn in glob(in_path + "*." + ext):
        try:
            metadata = parse_metadata(os.path.basename(fn))
            plate_prefix = metadata["plate_prefix"]
            plate = metadata["plate"]
            if plate not in dict_plate_prefixes.keys():
                dict_plate_prefixes[plate] = plate_prefix
            plates.append(plate)
            channels.append(f"A{metadata['A']}_C{metadata['C']}")
        except IndexError:
            print("IndexError for ", fn)
            pass
    plates = sorted(list(set(plates)))
    channels = sorted(list(set(channels)))
    print("Find all plates/channels in", in_path + "*." + ext)
    print(f"Plates:   {plates}")
    print(f"Channels: {channels}")

    # FIXME: remove hard-coded default (for None)
    if path_dict_channels is None:
        dict_channels = {
            "A01_C01": dict(label="DAPI", colormap="00FFFF", start=110),
            "A01_C02": dict(label="nanog", colormap="FF00FF", start=115),
            "A02_C03": dict(label="Lamin B1", colormap="FFFF00", start=115),
        }
    else:
        try:
            with open(path_dict_channels, "r") as json_file:
                dict_channels = json.load(json_file)
        except FileNotFoundError:
            raise Exception(
                "ERROR in create_zarr_structure: "
                f"{path_dict_channels} missing."
            )
        except TypeError:
            raise Exception(
                "ERROR in create_zarr_structure: "
                f"{path_dict_channels} has wrong type "
                "(probably a None instead of a string)."
            )

    # Check that all channels are in the allowed_channels
    if not set(channels).issubset(set(dict_channels.keys())):
        msg = "ERROR in create_zarr_structure\n"
        msg += f"channels: {channels}\n"
        msg += f"allowed_channels: {dict_channels.keys()}\n"
        raise Exception(msg)

    # Sort channels according to allowed_channels, and assign increasing index
    # actual_channels is a list of entries like A01_C01"
    actual_channels = []
    print("-" * 80)
    print("actual_channels")
    print("-" * 80)
    for ind_ch, ch in enumerate(channels):
        actual_channels.append([ch])
        print(actual_channels[-1])
    print("-" * 80)

    well = []

    zarrurls = {"plate": [], "well": []}

    for plate in plates:

        # Define plate zarr
        group_plate = zarr.group(out_path + f"{plate}.zarr")
        zarrurls["plate"].append(out_path + f"{plate}.zarr")

        # Identify all wells
        plate_prefix = dict_plate_prefixes[plate]
        wells = [
            parse_metadata(os.path.basename(fn))["well"]
            for fn in glob(f"{in_path}{plate_prefix}_*.{ext}")
        ]
        wells = sorted(list(set(wells)))

        # Verify that all wells have all channels
        for well in wells:
            well_channels = []
            glob_string = f"{in_path}{plate_prefix}_{well}*.{ext}"
            for fn in glob(glob_string):
                try:
                    metadata = parse_metadata(os.path.basename(fn))
                    well_channels.append(f"A{metadata['A']}_C{metadata['C']}")
                except IndexError:
                    print(f"Skipping {fn}")
            well_channels = sorted(list(set(well_channels)))
            if well_channels != actual_channels:
                raise Exception(
                    f"ERROR: well {well} in plate {plate} (prefix: "
                    f"{plate_prefix}) has missing channels.\n"
                    f"Expected: {actual_channels}\n"
                    f"Found: {well_channels}.\n"
                    f"[glob_string: {glob_string}]"
                )

        well_rows_columns = [
            ind for ind in sorted([(n[0], n[1:]) for n in wells])
        ]

        group_plate.attrs["plate"] = {
            "acquisitions": [
                # FIXME this should not be within "for plate in plates"!
                # {"id": 1, "name": plate}  # new version
                {"id": id_, "name": name}
                for id_, name in enumerate(plates)
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
                "id": 1,  # FIXME does this depend on the plate number?
                "name": "TBD",
                "version": "0.4",
                "channels": [
                    {
                        # "active": true,  # how to write it in python?
                        "coefficient": 1,
                        "color": dict_channels[channel]["colormap"],
                        "family": "linear",
                        # "inverted": false, # how to write it in python?
                        "label": channel[2],
                        "window": {
                            "min": 0,
                            "max": 65535,
                            "start": dict_channels[channel]["start"],
                            "end": dict_channels[channel]["end"],
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
