---
hide:
  - toc
---

# Welcome to Fractal Command-line Client's documentation!

Fractal is a framework to process high content imaging data at scale and prepare it for interactive visualization.

> This project is under active development 🔨. If you need help or found a bug, **open an issue [here](https://github.com/fractal-analytics-platform/fractal/issues/new)**.

Fractal provides distributed workflows that convert TBs of image data into [OME-Zar](https://ngff.openmicroscopy.org) files.
The platform then processes the 3D image data by applying tasks like illumination correction, maximum intensity projection, 3D segmentation using [cellpose](https://cellpose.readthedocs.io) and measurements using [napari workflows](https://github.com/haesleinhuepf/napari-workflows).
The pyramidal OME-Zarr files enable interactive visualization in the napari viewer.

This documentation concerns the **Fractal Command-line Client**. Find more information about Fractal in general and the other repositories at the [Fractal home page](https://fractal-analytics-platform.github.io).
