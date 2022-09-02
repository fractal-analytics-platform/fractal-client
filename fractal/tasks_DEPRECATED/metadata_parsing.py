"""
Copyright 2022 (C) Friedrich Miescher Institute for Biomedical Research and
University of Zurich

Original authors:
Joel Lüthi  <joel.luethi@fmi.ch>

This file is part of Fractal and was originally developed by eXact lab S.r.l.
<exact-lab.it> under contract with Liberali Lab from the Friedrich Miescher
Institute for Biomedical Research and Pelkmans Lab from the University of
Zurich.
"""
import warnings

import numpy as np
import pandas as pd
from defusedxml import ElementTree


def parse_yokogawa_metadata(mrf_path, mlf_path):
    """
    Parse Yokogawa CV7000 metadata files and prepare site-level metadata

    :param mrf_path: Full path to MeasurementDetail.mrf metadata file
    :type mrf_path: Union(str, pathlib.Path)
    :param mlf_path: Full path to MeasurementData.mlf metadata file
    :type mlf_path: Union(str, pathlib.Path)
    """
    mrf_frame, mlf_frame, error_count = read_metadata_files(mrf_path, mlf_path)

    per_site_parameters = [
        "x_micrometer",
        "y_micrometer",
        "pixel_size_x",
        "pixel_size_y",
        "x_pixel",
        "y_pixel",
        "bit_depth",
    ]
    grouping_params = ["well_id", "field_id"]
    grouped_sites = mlf_frame.loc[
        :, grouping_params + per_site_parameters
    ].groupby(by=grouping_params)
    check_grouped_sites_consistency(grouped_sites, per_site_parameters)
    site_metadata = grouped_sites.mean()

    # Cast image pixel sizes & bit depth to int
    site_metadata = site_metadata.astype(
        {
            "x_pixel": "int",
            "y_pixel": "int",
            "bit_depth": "int",
        }
    )

    # Absolute Z positions are not saved by the Yokogawa,
    # only relative positions to the autofocus
    site_metadata["z_micrometer"] = 0

    site_metadata = pd.concat(
        [
            site_metadata,
            get_z_steps(mlf_frame),
            get_earliest_time_per_site(mlf_frame),
        ],
        axis=1,
    )

    if error_count > 0:
        print(
            f"Succesfully parsed {len(site_metadata)} sites, could not "
            f"parse {error_count} sites due to errors (see warnings)."
        )
    total_files = len(mlf_frame)
    # TODO: Check whether the total_files correspond to the number of
    # relevant input images in the input folder. Returning it for now
    # Maybe return it here for further checks and produce a warning if it does
    # not match

    return site_metadata, total_files


def read_metadata_files(mrf_path, mlf_path):
    # parsing of mrf & mlf files are based on the
    # yokogawa_image_collection_task v0.5 in drogon, written by Dario Vischi.
    # https://github.com/fmi-basel/job-system-workflows/blob/00bbf34448972d27f258a2c28245dd96180e8229/src/gliberal_workflows/tasks/yokogawa_image_collection_task/versions/version_0_5.py  # noqa
    # Now modified for Fractal use

    mrf_frame = read_mrf_file(mrf_path)
    # TODO: filter_position & filter_wheel_position are parsed, but not
    # processed further. Figure out how to save them as relevant metadata for
    # use e.g. during illumination correction

    mlf_frame, error_count = read_mlf_file(mlf_path, mrf_frame)
    # TODO: Time points are parsed as part of the mlf_frame, but currently not
    # processed further. Once we tackle time-resolved data, parse from here.

    return mrf_frame, mlf_frame, error_count


def read_mrf_file(mrf_path):

    # Prepare mrf dataframe
    mrf_columns = [
        "channel_id",
        "horiz_pixel_dim",
        "vert_pixel_dim",
        "camera_no",
        "bit_depth",
        "horiz_pixels",
        "vert_pixels",
        "filter_wheel_position",
        "filter_position",
        "shading_corr_src",
    ]
    mrf_frame = pd.DataFrame(columns=mrf_columns)

    mrf_xml = ElementTree.parse(mrf_path).getroot()
    # Read mrf file
    ns = {"bts": "http://www.yokogawa.co.jp/BTS/BTSSchema/1.0"}
    for channel in mrf_xml.findall("bts:MeasurementChannel", namespaces=ns):
        mrf_frame.loc[channel.get("{%s}Ch" % ns["bts"])] = [
            channel.get("{%s}Ch" % ns["bts"]),
            float(channel.get("{%s}HorizontalPixelDimension" % ns["bts"])),
            float(channel.get("{%s}VerticalPixelDimension" % ns["bts"])),
            int(channel.get("{%s}CameraNumber" % ns["bts"])),
            int(channel.get("{%s}InputBitDepth" % ns["bts"])),
            int(channel.get("{%s}HorizontalPixels" % ns["bts"])),
            int(channel.get("{%s}VerticalPixels" % ns["bts"])),
            int(channel.get("{%s}FilterWheelPosition" % ns["bts"])),
            int(channel.get("{%s}FilterPosition" % ns["bts"])),
            channel.get("{%s}ShadingCorrectionSource" % ns["bts"]),
        ]

    return mrf_frame


def read_mlf_file(mlf_path, mrf_frame):
    # Docstring TBD
    mlf_xml = ElementTree.parse(mlf_path).getroot()
    ns = {"bts": "http://www.yokogawa.co.jp/BTS/BTSSchema/1.0"}

    # Measure number of lines
    def blocks(fh, size=65536):
        while True:
            block = fh.read(size)
            if not block:
                break
            yield block

    mlf_entries = mlf_xml.findall("bts:MeasurementRecord", namespaces=ns)
    nb_lines = len(mlf_entries)

    # Prepare mlf dataframe
    mlf_columns = [
        "type",
        "well_id",
        "column",
        "row",
        "time_point",
        "field_id",
        "z_index",
        "timeline_id",
        "action_id",
        "action",
        "x_micrometer",
        "y_micrometer",
        "z_micrometer",
        "x_pixel",
        "y_pixel",
        "pixel_size_x",
        "pixel_size_y",
        "bit_depth",
        "width",
        "height",
        "channel_id",
        "camera_no",
        "file_name",
        "time",
    ]
    mlf_frame = pd.DataFrame(columns=mlf_columns, index=range(0, nb_lines))

    mrf_channel_indices = {
        row.channel_id: idx
        for idx, (_, row) in enumerate(mrf_frame.iterrows())
    }

    error_count = 0
    for idx, record in enumerate(mlf_entries):
        rec_type = record.get("{%s}Type" % ns["bts"])

        if rec_type == "ERR":
            warnings.warn(
                "When parsing the yokogawa metadata, "
                f"found an 'ERR' entry at line '{idx+3}' in the MLF "
                f"meta file! The error was: '{record.text}'."
                "The entry is skipped."
            )
            error_count += 1
            continue
        elif rec_type != "IMG":
            warnings.warn(
                "When parsing the yokogawa metadata, "
                f"Found unexpected '{rec_type}' entry at line '{idx+3}'"
                "in the MLF meta file! The entry is skipped."
            )
            error_count += 1
            continue

        channel_id = record.get("{%s}Ch" % ns["bts"])
        x_micrometer = float(record.get("{%s}X" % ns["bts"]))
        # we mirror the y coordinate to fit with the field layout
        y_micrometer = -float(record.get("{%s}Y" % ns["bts"]))
        z_micrometer = float(record.get("{%s}Z" % ns["bts"]))

        well_row_id = record.get("{%s}Row" % ns["bts"])
        well_col_id = record.get("{%s}Column" % ns["bts"])
        well_id = chr(64 + int(well_row_id)) + str(well_col_id).zfill(2)
        # Convert all times to UTC time zone to avoid later timezone handling
        time = pd.to_datetime(record.get("{%s}Time" % ns["bts"]), utc=True)

        bit_depth = np.nan
        width = np.nan
        height = np.nan
        camera_no = np.nan
        pixel_size_x = np.nan
        pixel_size_y = np.nan
        if rec_type == "IMG":
            mrf_idx = mrf_channel_indices[channel_id]
            pixel_size_x = mrf_frame.iat[mrf_idx, 1]
            pixel_size_y = mrf_frame.iat[mrf_idx, 2]
            bit_depth = int(mrf_frame.iat[mrf_idx, 4])
            width = int(mrf_frame.iat[mrf_idx, 5])
            height = int(mrf_frame.iat[mrf_idx, 6])
            camera_no = int(mrf_frame.iat[mrf_idx, 3])

        # we use iat[] here for (significant) performance reasons
        mlf_frame.iat[idx, 0] = rec_type
        mlf_frame.iat[idx, 1] = well_id
        mlf_frame.iat[idx, 2] = int(well_col_id)
        mlf_frame.iat[idx, 3] = int(well_row_id)
        mlf_frame.iat[idx, 4] = int(record.get("{%s}TimePoint" % ns["bts"]))
        mlf_frame.iat[idx, 5] = int(record.get("{%s}FieldIndex" % ns["bts"]))
        mlf_frame.iat[idx, 6] = int(record.get("{%s}ZIndex" % ns["bts"]))
        mlf_frame.iat[idx, 7] = int(
            record.get("{%s}TimelineIndex" % ns["bts"])
        )
        mlf_frame.iat[idx, 8] = int(record.get("{%s}ActionIndex" % ns["bts"]))
        mlf_frame.iat[idx, 9] = record.get("{%s}Action" % ns["bts"])
        mlf_frame.iat[idx, 10] = x_micrometer
        mlf_frame.iat[idx, 11] = y_micrometer
        mlf_frame.iat[idx, 12] = z_micrometer
        mlf_frame.iat[idx, 13] = width
        mlf_frame.iat[idx, 14] = height
        mlf_frame.iat[idx, 15] = pixel_size_x
        mlf_frame.iat[idx, 16] = pixel_size_y

        mlf_frame.iat[idx, 17] = bit_depth
        mlf_frame.iat[idx, 18] = width
        mlf_frame.iat[idx, 19] = height
        mlf_frame.iat[idx, 20] = channel_id
        mlf_frame.iat[idx, 21] = camera_no
        mlf_frame.iat[idx, 22] = record.text  # file_name
        mlf_frame.iat[idx, 23] = time  # file_name

    mlf_frame = mlf_frame.dropna(thresh=(len(mlf_frame.columns)))
    return mlf_frame, error_count


def calculate_steps(site_series: pd.Series):
    # site_series is the z_micrometer series for a given site of a given
    # channel. This function calculates the step size in Z

    # First diff is always NaN because there is nothing to compare it to
    steps = site_series.diff().dropna()
    if not steps.std().sum() == 0.0:
        raise Exception(
            "When parsing the Yokogawa mlf file, some sites "
            "had varying step size in Z. "
            "That is not supported for the OME-Zarr parsing"
        )
    return steps.mean()


def get_z_steps(mlf_frame):
    # Process mlf_frame to extract Z information (pixel size & steps).
    # Run checks on consistencies & return site-based z step dataframe
    # Group by well, field & channel
    grouped_sites_z = (
        mlf_frame.loc[
            :,
            ["well_id", "field_id", "action_id", "channel_id", "z_micrometer"],
        ]
        .set_index(["well_id", "field_id", "action_id", "channel_id"])
        .groupby(level=[0, 1, 2, 3])
    )

    # If there is only 1 Z step, set the Z spacing to the count of planes => 1
    if grouped_sites_z.count()["z_micrometer"].max() == 1:
        z_data = grouped_sites_z.count().groupby(["well_id", "field_id"])
    else:
        # Group the whole site (combine channels), because Z steps need to be
        # consistent between channels for OME-Zarr.
        z_data = grouped_sites_z.apply(calculate_steps).groupby(
            ["well_id", "field_id"]
        )
    if not z_data.std().sum().sum() == 0.0:
        raise Exception(
            "When parsing the Yokogawa mlf file, channels had "
            "varying step size in Z. "
            "That is not supported for the OME-Zarr parsing"
        )

    # Ensure that channels have the same number of z planes and
    # reduce it to one value.
    # Only check if there is more than one channel available
    if any(
        grouped_sites_z.count().groupby(["well_id", "field_id"]).count() > 1
    ):
        if any(
            grouped_sites_z.count()
            .groupby(["well_id", "field_id"])
            .std()
            .sum()
            != 0
        ):
            raise Exception(
                "When parsing the Yokogawa mlf file, channels had "
                "varying number of z planes."
                "That is not supported for the OME-Zarr parsing"
            )
    z_steps = (
        grouped_sites_z.count()
        .groupby(["well_id", "field_id"])
        .mean()
        .astype(int)
    )

    # Combine the two dataframes
    z_frame = pd.concat([z_data.mean(), z_steps], axis=1)
    z_frame.columns = ["pixel_size_z", "z_pixel"]
    return z_frame


def check_grouped_sites_consistency(grouped_sites, per_site_parameters):
    # Check that stage positions don't vary within a given site
    # Same for pixel sizes
    # Only relevant when a site has multiple entries
    if grouped_sites.count().min().min() > 1:
        if not grouped_sites.std().sum().sum() == 0.0:
            raise Exception(
                ""
                "When parsing the Yokogawa MeasurementData.mlf file, "
                f"some of the parameters {per_site_parameters} varied within "
                "the site. That is not supported for the OME-Zarr parsing"
            )


def get_earliest_time_per_site(mlf_frame) -> pd.DataFrame:
    # Get the time information per site
    # Because a site will contain time information for each plane
    # of each channel, we just return the earliest time infromation
    # per site.
    return mlf_frame.groupby(["well_id", "field_id"]).min()["time"]
