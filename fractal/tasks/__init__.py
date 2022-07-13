__FRACTAL_MANIFEST__ = [
    {
        "resource_type": "core task",
        "name": "dummy",
        "module": f"{__name__}.dummy:dummy",
        "input_type": "Any",
        "output_type": "None",
    },
    {
        "resource_type": "core task",
        "name": "Create OME-ZARR structure",
        "module": f"{__name__}.create_zarr_structure:create_zarr_structure",
        "input_type": "image",
        "output_type": "zarr",
        "default_parameters": {
            "num_levels": 2,
            "channel_parameters": {
                "A01_C01": {
                    "label": "DAPI",
                    "colormap": "00FFFF",
                    "start": 0,
                    "end": 700,
                },
                "A01_C02": {
                    "label": "nanog",
                    "colormap": "FF00FF",
                    "start": 0,
                    "end": 180,
                },
                "A02_C03": {
                    "label": "Lamin B1",
                    "colormap": "FFFF00",
                    "start": 0,
                    "end": 1500,
                },
            },
        },
    },
]
