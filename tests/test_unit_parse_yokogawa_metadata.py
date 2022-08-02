import os

from fractal.tasks.metadata_parsing import parse_yokogawa_metadata


def test_parse_yokogawa_metadata():
    testdir = os.path.dirname(__file__)
    path = f"{testdir}/data/png/"
    metadata, num = parse_yokogawa_metadata(
        path + "MeasurementDetail.mrf", path + "MeasurementData.mlf"
    )
    assert num == 4


if __name__ == "__main__":
    test_parse_yokogawa_metadata()
