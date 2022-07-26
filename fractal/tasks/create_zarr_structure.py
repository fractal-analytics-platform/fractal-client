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
import os
from glob import glob
from pathlib import Path
from typing import Any
from typing import Dict
from typing import Iterable
from typing import Optional

import pandas as pd
import zarr
from anndata.experimental import write_elem

from fractal.tasks.lib_parse_filename_metadata import parse_metadata
from fractal.tasks.lib_regions_of_interest import prepare_FOV_ROI_table
from fractal.tasks.metadata_parsing import parse_yokogawa_metadata


def create_zarr_structure(
    *,
    input_paths: Iterable[Path],
    output_path: Path,
    channel_parameters: Dict[str, Any],
    num_levels: int = 2,
    coarsening_xy: int = 2,
    metadata_table: str = "mrf_mlf",
    metadata: Optional[Dict[str, Any]] = None,
):

    """
    Create (and store) the zarr folder, without reading or writing data.

    1. Find plates
        For each folder in input paths:
        * glob image files
        * parse metadata from image filename to identify plates
        * identify populated channels

    2. Create a ZARR for each plate
        For each plate:
        * parse mlf metadata
        * identify wells and field of view (FOV)
        * create FOV ZARR
        * verify that channels are uniform (i.e., same channels)

    :param in_paths: list of image directories
    :type in_path: list
    :param out_path: path for output zarr files
    :type out_path: str
    :param ext: extension of images (e.g. tiff, png, ..)
    :param path_dict_channels: FIXME
    :type path_dict_channels: str
    :param num_levels: number of coarsening levels in the pyramid
    :type num_levels: int
    FIXME
    """

    # Preliminary checks on metadata_table
    if metadata_table != "mrf_mlf" and not isinstance(
        metadata_table, pd.core.frame.DataFrame
    ):
        raise Exception(
            "ERROR: metadata_table must be a known string or a "
            "pandas DataFrame}"
        )
    if metadata_table != "mrf_mlf":
        raise NotImplementedError(
            "We currently only support "
            'metadata_table="mrf_mlf", '
            f"and not {metadata_table}"
        )

    # Identify all plates and all channels, across all input folders
    plates = []
    channels = None
    dict_plate_paths = {}
    dict_plate_prefixes = {}

    # FIXME
    # find a smart way to remove it
    ext_glob_pattern = input_paths[0].name

    for in_path in input_paths:
        input_filename_iter = in_path.parent.glob(in_path.name)

        tmp_channels = []
        tmp_plates = []
        for fn in input_filename_iter:
            try:
                metadata = parse_metadata(fn.name)
                plate_prefix = metadata["plate_prefix"]
                plate = metadata["plate"]
                if plate not in dict_plate_prefixes.keys():
                    dict_plate_prefixes[plate] = plate_prefix
                tmp_plates.append(plate)
                tmp_channels.append(f"A{metadata['A']}_C{metadata['C']}")
            except IndexError:
                print("IndexError for ", fn)
                pass
        tmp_plates = sorted(list(set(tmp_plates)))
        tmp_channels = sorted(list(set(tmp_channels)))

        info = (
            f"Listing all plates/channels from {in_path.as_posix()}\n"
            f"Plates:   {tmp_plates}\n"
            f"Channels: {tmp_channels}\n"
        )

        # Check that only one plate is found
        if len(tmp_plates) > 1:
            raise Exception(f"{info}ERROR: {len(tmp_plates)} plates detected")
        plate = tmp_plates[0]

        # If plate already exists in other folder, add suffix
        if plate in plates:
            ind = 1
            new_plate = f"{plate}_{ind}"
            while new_plate in plates:
                new_plate = f"{plate}_{ind}"
                ind += 1
            print(
                f"WARNING: {plate} already exists, renaming it as {new_plate}"
            )
            plates.append(new_plate)
            dict_plate_prefixes[new_plate] = dict_plate_prefixes[plate]
            plate = new_plate
        else:
            plates.append(plate)

        # Check that channels are the same as in previous plates
        if channels is None:
            channels = tmp_channels[:]
        else:
            if channels != tmp_channels:
                raise Exception(
                    f"ERROR\n{info}\nERROR: expected channels " "{channels}"
                )

        # Update dict_plate_paths
        dict_plate_paths[plate] = in_path.parent

    # Check that all channels are in the allowed_channels
    if not set(channels).issubset(set(channel_parameters.keys())):
        msg = "ERROR in create_zarr_structure\n"
        msg += f"channels: {channels}\n"
        msg += f"allowed_channels: {channel_parameters.keys()}\n"
        raise Exception(msg)

    # Sort channels according to allowed_channels, and assign increasing index
    # actual_channels is a list of entries like A01_C01"
    actual_channels = []
    for ind_ch, ch in enumerate(channels):
        actual_channels.append(ch)
    print(f"actual_channels: {actual_channels}")

    zarrurls = {"plate": [], "well": []}

    ################################################################
    for plate in plates:
        # Define plate zarr
        zarrurl = f"{output_path.as_posix()}/{plate}.zarr"
        print(f"Creating {zarrurl}")
        group_plate = zarr.group(zarrurl)
        zarrurls["plate"].append(zarrurl)

        # Obtain FOV-metadata dataframe
        try:
            # FIXME
            # Find a smart way to include these metadata files in the dataset
            # e.g., as resources
            if metadata_table == "mrf_mlf":
                mrf_path = f"{in_path}MeasurementDetail.mrf"
                mlf_path = f"{in_path}MeasurementData.mlf"
                site_metadata, total_files = parse_yokogawa_metadata(
                    mrf_path, mlf_path
                )
                has_mrf_mlf_metadata = True

                # Extract pixel sizes
                pixel_size_z = site_metadata["pixel_size_z"][0]
                pixel_size_y = site_metadata["pixel_size_y"][0]
                pixel_size_x = site_metadata["pixel_size_x"][0]
                # Extract bit_depth #FIXME
                # bit_depth = site_metadata["bit_depth"][0]
                # if bit_depth == 8:
                #    dtype
        except FileNotFoundError:
            print("Missing metadata files")
            has_mrf_mlf_metadata = False

            pixel_size_x = pixel_size_y = pixel_size_z = 1

        # Identify all wells
        plate_prefix = dict_plate_prefixes[plate]
        in_path = dict_plate_paths[plate]

        plate_image_iter = glob(f"{in_path}/{plate_prefix}_{ext_glob_pattern}")

        wells = [
            parse_metadata(os.path.basename(fn))["well"]
            for fn in plate_image_iter
        ]
        wells = sorted(list(set(wells)))

        # Verify that all wells have all channels
        for well in wells:
            well_image_iter = glob(
                f"{in_path}/{plate_prefix}_{well}{ext_glob_pattern}"
            )
            well_channels = []
            for fn in well_image_iter:
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
                )

        well_rows_columns = [
            ind for ind in sorted([(n[0], n[1:]) for n in wells])
        ]

        group_plate.attrs["plate"] = {
            "acquisitions": [{"id": 1, "name": plate}],
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
                "images": [{"path": "0"}],
                "version": "0.3",
            }

            group_FOV = group_well.create_group("0/")  # noqa: F841
            zarrurls["well"].append(
                output_path.as_posix() + f"/{plate}.zarr/{row}/{column}/0/"
            )

            group_FOV.attrs["multiscales"] = [
                {
                    "version": "0.3",
                    "axes": [
                        {"name": "c", "type": "channel"},
                        {
                            "name": "z",
                            "type": "space",
                            "unit": "micrometer",
                        },
                        {
                            "name": "y",
                            "type": "space",
                            "unit": "micrometer",
                        },
                        {
                            "name": "x",
                            "type": "space",
                            "unit": "micrometer",
                        },
                    ],
                    "datasets": [
                        {
                            "path": f"{ind_level}",
                            "coordinateTransformations": [
                                {
                                    "type": "scale",
                                    "scale": [
                                        pixel_size_z,
                                        pixel_size_y
                                        * coarsening_xy**ind_level,
                                        pixel_size_x
                                        * coarsening_xy**ind_level,
                                    ],
                                }
                            ],
                        }
                        for ind_level in range(num_levels)
                    ],
                }
            ]

            group_FOV.attrs["omero"] = {
                "id": 1,  # FIXME does this depend on the plate number?
                "name": "TBD",
                "version": "0.4",
                "channels": [
                    {
                        "active": True,
                        "coefficient": 1,
                        "color": channel_parameters[channel]["colormap"],
                        "family": "linear",
                        "inverted": False,
                        "label": channel_parameters[channel]["label"],
                        "window": {
                            "min": 0,
                            "max": 65535,
                            "start": channel_parameters[channel]["start"],
                            "end": channel_parameters[channel]["end"],
                        },
                    }
                    for channel in actual_channels
                ],
            }

            # FIXME
            if has_mrf_mlf_metadata:
                # Prepare and write anndata table of FOV ROIs
                FOV_ROIs_table = prepare_FOV_ROI_table(
                    site_metadata.loc[f"{row+column}"],
                )
                group_tables = group_FOV.create_group("tables/")  # noqa: F841
                write_elem(group_tables, "FOV_ROI_table", FOV_ROIs_table)

    return zarrurls, actual_channels


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(prog="create_zarr_structure")
    parser.add_argument(
        "-i",
        "--in_paths",
        help="list of directories containing the input files",
        nargs="+",
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
    parser.add_argument(
        "-cxy",
        "--coarsening_xy",
        type=int,
        help="FIXME",
    )
    parser.add_argument(
        "-c",
        "--path_dict_channels",
        type=str,
        help="path of channel dictionary",
    )
    args = parser.parse_args()
    create_zarr_structure(
        in_paths=args.in_paths,
        out_path=args.out_path,
        ext=args.ext,
        num_levels=args.num_levels,
        coarsening_xy=args.coarsening_xy,
        path_dict_channels=args.path_dict_channels,
        # metadata_table=args.metadata_table,   #FIXME
    )
