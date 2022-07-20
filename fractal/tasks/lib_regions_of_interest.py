import anndata as ad
import numpy as np


def prepare_ROIs_table(df, image_size=None):
    if image_size is None:
        raise Exception("Missing image_size arg in prepare_ROIs_table")

    df["x_micrometer"] -= df["x_micrometer"].min()
    df["y_micrometer"] -= df["y_micrometer"].min()
    df["z_micrometer"] -= df["z_micrometer"].min()

    df["len_x_micrometer"] = image_size["x"] * df["pixel_size_x"]
    df["len_y_micrometer"] = image_size["y"] * df["pixel_size_y"]
    df["len_z_micrometer"] = df["pixel_size_z"]

    df.drop("bit_depth", inplace=True, axis=1)

    df = df.astype(np.float32)

    adata = ad.AnnData(X=df, dtype=np.float32)
    adata.obs_names = [f"FOV_{i+1:d}" for i in range(len(df.index))]
    adata.var_names = df.columns

    return adata
