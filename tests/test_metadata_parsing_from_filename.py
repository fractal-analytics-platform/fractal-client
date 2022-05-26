from fractal.tasks.create_zarr_structure import metadata

f1 = "20200812-CardiomyocyteDifferentiation14-Cycle1"
f1 += "_B03_T0001F036L01A01Z18C01.png"
f2 = "210305NAR005AAN_210416_164828_B11_T0001F006L01A04Z14C01.tif"
f3 = "220304_172545_220304_175557_L06_T0277F004L277A04Z07C04.tif"

p1 = "20200812-CardiomyocyteDifferentiation14-Cycle1"
p2 = "210305NAR005AAN"
p3 = "RS220304172545"


def test_metadata():
    assert metadata(f1)["plate"] == p1
    assert metadata(f2)["plate"] == p2
    assert metadata(f3)["plate"] == p3
