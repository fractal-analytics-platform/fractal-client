import os

import pytest

from fractal.tasks.metadata_parsing import parse_yokogawa_metadata


@pytest.mark.skip("Skip due to issue 130")
def test_parse_yokogawa_metadata():
    # General variables and paths (relative to mwe_fractal folder)
    rootdir, relative_dir = os.getcwd().split("mwe_fractal")
    if relative_dir != "":
        raise Exception(
            "ERROR: this test has hard-coded paths, it should be "
            "run from mwe_fractal folder"
        )

    path = f"{rootdir}mwe_fractal/tests/data/png/"
    metadata, num = parse_yokogawa_metadata(
        path + "MeasurementDetail.mrf", path + "MeasurementData.mlf"
    )
    assert num == 2


if __name__ == "__main__":
    test_parse_yokogawa_metadata()
