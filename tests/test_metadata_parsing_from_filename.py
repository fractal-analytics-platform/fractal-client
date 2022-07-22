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
import pytest

from fractal.tasks.lib_parse_filename_metadata import parse_metadata

f1 = "20200812-CardiomyocyteDifferentiation14-Cycle1"
f1 += "_B03_T0001F036L01A01Z18C01.png"
f2 = "210305NAR005AAN_210416_164828_B11_T0001F006L01A04Z14C01.tif"
f3 = "220304_172545_220304_175557_L06_T0277F004L277A04Z07C04.tif"
f4 = "220517CS001XXXIl_220715_151525_D05_T0001F001L01A01Z01C04.tif"

p1 = "20200812-CardiomyocyteDifferentiation14-Cycle1"
p2 = "210305NAR005AAN"
p3 = "RS220304172545"
p4 = "220517CS001XXXIl"

A1 = "01"
A2 = "04"
A3 = "04"
A4 = "01"

C1 = "01"
C2 = "01"
C3 = "04"
C4 = "04"

parameters = [
    (f1, p1, A1, C1),
    (f2, p2, A2, C2),
    (f3, p3, A3, C3),
    (f4, p4, A4, C4),
]


@pytest.mark.parametrize("filename,plate,A,C", parameters)
def test_metadata(filename, plate, A, C):
    metadata = parse_metadata(filename)
    assert metadata["plate"] == plate
    assert metadata["A"] == A
    assert metadata["C"] == C
