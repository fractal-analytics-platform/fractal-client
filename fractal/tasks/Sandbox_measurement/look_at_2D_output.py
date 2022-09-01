import anndata as ad
from devtools import debug


table_path = (
    "/data/active/fractal/tests/"
    "Temporary_data_UZH_1_well_2x2_sites/"
    "20200812-CardiomyocyteDifferentiation14-Cycle1_mip.zarr"
    "/B/03/0/tables/nuclei"
)
adata = ad.read_zarr(table_path)
num_labels = len(list(adata.obs_names))
num_labels_unique = len(set(list(adata.obs_names)))

debug(adata)
print()
debug(adata.var)
print()
debug(type(adata.obs))
print()
debug(adata.obs.dtypes)
print()
debug(num_labels)
print()

if not num_labels == num_labels_unique:
    raise Exception("Non-unique labels")
