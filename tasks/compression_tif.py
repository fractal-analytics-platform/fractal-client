import glob
import json
import os
import sys

import tifffile
from PIL import Image


def metadata(filename):

    with tifffile.TiffFile(filename) as tif:
        tif_tags = {}
        for tag in tif.pages[0].tags.values():
            name, value = tag.name, tag.value
            tif_tags[name] = value
    return tif_tags


def compress(in_path, out_path, delete_in):

    if not os.path.exists(out_path):
        os.makedirs(out_path)

    img_metadata = {}
    file_list = []

    if isinstance(in_path, str):
        os.chdir(in_path)
        for filename in glob.glob("*.tif"):
            with Image.open(filename) as image:
                image.save(
                    os.path.join(out_path, filename.strip(".tif")+"_compressed.tif"),
                    compression="tiff_lzw",
                )
            file_list.append(filename)
            img_metadata[filename] = metadata(filename)

        if delete_in=="True":
            for f in file_list:
                try:
                    os.remove(f)
                except OSError as e:
                    print("Error: %s : %s" % (f, e.strerror))
        


        with open(out_path + "/" + "metadata.json", "w") as f:
            json.dump(img_metadata, f)


if __name__ == "__main__":
    in_path = sys.argv[1]
    out_path = sys.argv[2]
    delete_in = sys.argv[3]
    compress(in_path, out_path, delete_in)
