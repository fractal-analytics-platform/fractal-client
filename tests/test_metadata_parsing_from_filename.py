"""
Copyright 2022 (C) Friedrich Miescher Institute for Biomedical Research and
University of Zurich

Original authors:
Tommaso Comparin <tommaso.comparin@exact-lab.it>

This file is part of Fractal and was originally developed by eXact lab S.r.l.
<exact-lab.it> under contract with Liberali Lab from the Friedrich Miescher
Institute for Biomedical Research and Pelkmans Lab from the University of
Zurich.
"""

from fractal.tasks.lib_parse_filename_metadata import parse_metadata

f1 = "20200812-CardiomyocyteDifferentiation14-Cycle1"
f1 += "_B03_T0001F036L01A01Z18C01.png"
f2 = "210305NAR005AAN_210416_164828_B11_T0001F006L01A04Z14C01.tif"
f3 = "220304_172545_220304_175557_L06_T0277F004L277A04Z07C04.tif"

p1 = "20200812-CardiomyocyteDifferentiation14-Cycle1"
p2 = "210305NAR005AAN"
p3 = "RS220304172545"


def test_metadata():
    assert parse_metadata(f1)["plate"] == p1
    assert parse_metadata(f2)["plate"] == p2
    assert parse_metadata(f3)["plate"] == p3
