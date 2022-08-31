__FRACTAL_MANIFEST__ = [
    {
        "resource_type": "core task",
        "name": "dummy",
        "module": f"{__name__}.dummy:dummy",
        "input_type": "Any",
        "output_type": "None",
        "default_args": {
            "message": "dummy default",
            "index": 0,
            "needs_gpu": False,
        },
    },
    {
        "resource_type": "core task",
        "name": "Create OME-ZARR structure",
        "module": f"{__name__}.create_zarr_structure:create_zarr_structure",
        "input_type": "image",
        "output_type": "zarr",
        "default_args": {
            "needs_gpu": False,
            "num_levels": 2,
            "coarsening_xy": 2,
            "metadata_table": "mrf_mlf",
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
    {
        "name": "Yokogawa to Zarr",
        "resource_type": "core task",
        "input_type": "zarr",
        "output_type": "zarr",
        "module": f"{__name__}.yokogawa_to_zarr:yokogawa_to_zarr",
        "default_args": {"needs_gpu": False},
    },
    {
        "name": "Replicate Zarr structure",
        "resource_type": "core task",
        "input_type": "zarr",
        "output_type": "zarr",
        "module": f"{__name__}.replicate_zarr_structure:replicate_zarr_structure",
        "default_args": {
            "needs_gpu": False,
            "project_to_2D": True,
            "suffix": "mip",
        },
    },
    {
        "name": "Maximum Intensity Projection",
        "resource_type": "core task",
        "input_type": "zarr",
        "output_type": "zarr",
        "module": f"{__name__}.maximum_intensity_projection:maximum_intensity_projection",
        "default_args": {"needs_gpu": False},
    },
    {
        "name": "Per-FOV image labeling",
        "resource_type": "core task",
        "input_type": "zarr",
        "output_type": "zarr",
        "module": f"{__name__}.image_labeling:image_labeling",
        "default_args": {
            "labeling_channel": "A01_C01",
            "needs_gpu": False,  # FIXME
        },
    },
    {
        "name": "Whole-well image labeling",
        "resource_type": "core task",
        "input_type": "zarr",
        "output_type": "zarr",
        "module": f"{__name__}.image_labeling_whole_well:image_labeling_whole_well",
        "default_args": {
            "labeling_channel": "A01_C01",
            "needs_gpu": False,  # FIXME
        },
    },

]
