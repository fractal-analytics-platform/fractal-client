import os
import re
from glob import glob

import zarr

# from devtools import debug


def metadata(filename):
    """
    Extract metadata by parsing image filename, return a parameter dictionary.
    Three kinds of filenames are supported:
    1) Filenames from UZH:
       20200812-Cardio[...]Cycle1_B03_T0001F036L01A01Z18C01.png
       with plate name 20200812-Cardio[...]Cycle1
    2) Filenames from FMI, with successful barcode reading:
       210305NAR005AAN_210416_164828_B11_T0001F006L01A04Z14C01.tif
       with plate name 210305NAR005AAN
    3) Filenames from FMI, with failed barcode reading:
       yymmdd_hhmmss_210416_164828_B11_T0001F006L01A04Z14C01.tif
       with plate name RS{yymmddhhmmss}

    :param filename: name of the image
    :type filename: str
    """
    f = filename.rsplit(".", 1)[0]

    well = re.findall(r"_(.*)_T", f)[0].split("_")[-1]
    tmp_plate = f.split(f"_{well}_")[0]

    fields = tmp_plate.split("_")

    if (
        len(fields) == 4
        and len(fields[0]) == 6
        and len(fields[1]) == 6
        and len(fields[2]) == 6
    ):
        # FMI (failed barcode reading)
        scan_date, scan_time, img_date, img_time = fields[:]
        plate = f"RS{scan_date + scan_time}"
    elif len(fields) == 3:
        # FMI (correct barcode reading)
        barcode, img_date, img_time = fields[:]
        if len(img_date) != 6 or len(img_time) != 6:
            raise Exception(
                f"Failure in metadata parsing of {tmp_plate}, with"
                " img_date={img_date} and img_time={img_time}"
            )
        plate = barcode
    elif len(fields) == 1:
        # UZH
        plate = fields[0]

    site = re.findall(r"F(.*)L", f)[0]
    chl = re.findall(r"[0-9]C(.*)", f)[0].split(".")[0].split("_")[0]
    t_ind = re.findall(r"T(.*)F", f)[0]
    z_ind = re.findall(r"Z(.*)C", f)[0]

    result = dict(
        plate=plate, well=well, t_ind=t_ind, z_ind=z_ind, chl=chl, site=site
    )
    return result


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

    # Find all plates
    plate = []
    if not in_path.endswith("/"):
        in_path += "/"
    for i in glob(in_path + "*." + ext):
        try:
            plate.append(metadata(os.path.basename(i))["plate"])
        except IndexError:
            print("IndexError for ", i)
            pass
    plate_unique = set(plate)
    print("Find all plates in", in_path + "*." + ext)
    print(f"Plates: {plate_unique}")

    well = []

    zarrurls = {"plate": [], "well": []}

    # loop over plate, each plate could have n wells
    # debug(plate_unique)
    for plate in plate_unique:
        group_plate = zarr.group(out_path + f"{plate}.zarr")
        zarrurls["plate"].append(out_path + f"{plate}.zarr")
        well = [
            metadata(os.path.basename(fn))["well"]
            for fn in glob(in_path + f"{plate}_*." + ext)
        ]
        well_unique = set(well)

        well_rows_columns = [
            ind for ind in sorted([(n[0], n[1:]) for n in well_unique])
        ]

        group_plate.attrs["plate"] = {
            "acquisitions": [
                {"id": id_, "name": name}
                for id_, name in enumerate(plate_unique)
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
        # debug(well_rows_columns)

        for row, column in well_rows_columns:

            group_well = group_plate.create_group(f"{row}/{column}/")

            chl = [
                metadata(os.path.basename(fn))["chl"]
                for fn in glob(in_path + f"{plate}*_{row+column}*." + ext)
            ]
            chl_unique = set(chl)
            # debug(chl_unique)

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

    return zarrurls, chl_unique


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
