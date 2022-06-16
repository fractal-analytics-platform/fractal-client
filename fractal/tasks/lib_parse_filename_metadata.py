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

    # Parse filename for additional fields
    # Example filename: (...)_B03_T0001F001L01A01Z06C01.png
    T = re.findall(r"T(.*)F", f)[0]
    F = re.findall(r"F(.*)L", f)[0]
    L = re.findall(r"L(.*)A", f)[0]
    A = re.findall(r"A(.*)Z", f)[0]
    Z = re.findall(r"Z(.*)C", f)[0]
    C = re.findall(r"[0-9]C(.*)", f)[0].split(".")[0]

    result = dict(plate=plate, well=well, T=T, F=F, L=L, A=A, Z=Z, C=C)
    return result
