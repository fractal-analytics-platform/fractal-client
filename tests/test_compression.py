from io import BytesIO

from PIL import Image

from src.tasks.compression_tif import compress

in_path = ""
out_path = ""
start = "0"
end = "1"
delete_in = "False"


def test_compress(mocker):
    def create_test_image(self):
        file = BytesIO()
        image = Image.new("RGBA", size=(20, 20), color=(155, 0, 0))
        image.save(file, "tiff")
        file.name = "test.tiff"
        file.seek(0)
        rtn = Image.frombytes("L", (10, 10), file.read())
        return rtn

    mocker.patch("PIL.Image.open", create_test_image)

    mocker.patch("os.makedirs", return_value="")

    mocker.patch("glob.glob", return_value=["test.tiff"])

    buf = BytesIO()
    mocker.patch("os.path.join", return_value=buf)

    res = compress(in_path, start, end, out_path, delete_in)

    assert res == "tiff_lzw"
