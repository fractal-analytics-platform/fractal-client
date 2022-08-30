from pathlib import Path

from fractal.tasks.measurement import measurement


input_paths = [
    Path("/data/active/fractal/tests/Temporary_data_UZH_1_well_2x2_sites")
]
component = "20200812-CardiomyocyteDifferentiation14-Cycle1.zarr/B/03/0"
labeling_channel = "A01_C01"
level = 0
metadata = dict(coarsening_xy=2, channel_list=["A01_C01"])
workflow_file = "regionprops_from_existing_labels.yaml"
table_name = "nuclei"
whole_well = False


measurement(
    input_paths=input_paths,
    output_path=input_paths[0],
    metadata=metadata,
    component=component,
    labeling_channel=labeling_channel,
    level=level,
    workflow_file=workflow_file,
    table_name=table_name,
    whole_well=whole_well,
)
