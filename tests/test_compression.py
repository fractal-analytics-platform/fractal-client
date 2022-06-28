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

from io import BytesIO

from PIL import Image

from fractal.tasks.compress_tif import compress_tif

in_path = ""
out_path = ""
delete_input = False


def test_compress(mocker):
    def create_test_image(self):
        f = BytesIO()
        image = Image.new("RGBA", size=(20, 20), color=(155, 0, 0))
        image.save(f, "tiff")
        f.name = "test.tiff"
        f.seek(0)
        rtn = Image.frombytes("L", (10, 10), f.read())
        return rtn

    mocker.patch("PIL.Image.open", create_test_image)

    mocker.patch("os.makedirs", return_value="")

    mocker.patch("glob.glob", return_value=["test.tiff"])

    buf = BytesIO()
    mocker.patch("os.path.join", return_value=buf)

    res = compress_tif(in_path, out_path, delete_input=delete_input)

    assert res == (1, 0)
