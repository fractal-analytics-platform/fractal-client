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

import glob
import os

from PIL import Image


def compress_tif(in_path, out_path, delete_input=False):

    """
    Compress tiff files

    :param in_path: directory containing the input files
    :type in_path: str
    :param out_path: directory containing the output files
    :type out_path: str
    :param delete_input: delete input files
    :type delete_input: bool
    """

    # Sanitize input/output paths
    if not in_path.endswith("/"):
        in_path += "/"
    if not out_path.endswith("/"):
        out_path += "/"

    # Create output path, if needed
    if not os.path.exists(out_path):
        os.makedirs(out_path)

    num_img_compressed = 0
    num_img_deleted = 0
    for filename in glob.glob(in_path + "*.tif"):
        newfilename = os.path.join(out_path, os.path.basename(filename))

        # Save compressed image
        with Image.open(filename) as image:
            image.save(newfilename, format="tiff", compression="tiff_lzw")
        print(f"Raw:        {filename}\nCompressed: {newfilename}")
        num_img_compressed += 1

        # Delete raw image, if needed
        if delete_input:
            try:
                os.remove(filename)
                print(f"Deleted:    {filename}")
                num_img_deleted += 1
            except OSError as e:
                print("ERROR: %s : %s" % (filename, e.strerror))

        print()

    return num_img_compressed, num_img_deleted


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(prog="compress_tif")
    parser.add_argument("in_path", help="directory containing the input files")
    parser.add_argument(
        "out_path", help="directory containing the output files"
    )
    parser.add_argument(
        "-d",
        "--delete_input",
        action="store_true",
        help="Delete input files",
    )

    args = parser.parse_args()

    compress_tif(
        args.in_path,
        args.out_path,
        delete_input=args.delete_input,
    )
