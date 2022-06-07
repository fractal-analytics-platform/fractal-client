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
