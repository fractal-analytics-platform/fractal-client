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
import re


def parse_metadata(filename):
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

    if "/" in filename:
        raise Exception(
            "ERROR: parse_metadata may fail when filename "
            f'includes "/". Please check that {filename} is '
            "correct."
        )
    f = filename.rsplit(".", 1)[0]
    well = re.findall(r"_(.*)_T", f)[0].split("_")[-1]
    plate_prefix = f.split(f"_{well}_")[0]
    fields = plate_prefix.split("_")
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
                f"Failure in metadata parsing of {plate_prefix}, with"
                " img_date={img_date} and img_time={img_time}"
            )
        plate = barcode
    elif len(fields) == 1:
        # UZH
        plate = fields[0]

    # Parse filename for additional fields
    # Example of f_without_prefix: B03_T0001F001L01A01Z06C01.png
    f_without_prefix = f.split(plate_prefix + "_")[1]
    T = re.findall(r"_T(.*)F", f_without_prefix)[0]
    F = re.findall(rf"_T{T}F(.*)L", f_without_prefix)[0]
    L = re.findall(rf"_T{T}F{F}L(.*)A", f_without_prefix)[0]
    A = re.findall(rf"_T{T}F{F}L{L}A(.*)Z", f_without_prefix)[0]
    Z = re.findall(rf"_T{T}F{F}L{L}A{A}Z(.*)C", f_without_prefix)[0]
    C = re.findall(rf"_T{T}F{F}L{L}A{A}Z{Z}C(.*)", f_without_prefix)[0]

    result = dict(
        plate=plate,
        plate_prefix=plate_prefix,
        well=well,
        T=T,
        F=F,
        L=L,
        A=A,
        Z=Z,
        C=C,
    )
    return result
