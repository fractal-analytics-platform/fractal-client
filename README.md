# Fractal Client

[![PyPI version](https://img.shields.io/pypi/v/fractal-client?color=gree)](https://pypi.org/project/fractal-client/)

Fractal is a framework to process high content imaging data at scale and prepare it for interactive visualization.

Fractal provides distributed workflows that convert TBs of image data into OME-Zarr files. The platform then processes the 3D image data by applying tasks like illumination correction, maximum intensity projection, 3D segmentation using [cellpose](https://cellpose.readthedocs.io/en/latest/) and measurements using [napari workflows](https://github.com/haesleinhuepf/napari-workflows). The pyramidal OME-Zarr files enable interactive visualization in the napari viewer.

![Fractal_Overview](https://user-images.githubusercontent.com/18033446/190978261-2e7b57e9-72c7-443e-9202-15d233f8416d.jpg)


This is the main Fractal repository that contains the **Fractal client** and some examples. The **Fractal core tasks** to parse images and process OME-Zarr files can be found [here](https://github.com/fractal-analytics-platform/fractal-tasks-core). The **Fractal server** can be found [here](https://github.com/fractal-analytics-platform/fractal-server).

Example input data for Fractal can be found here: [10.5281/zenodo.7057076](https://doi.org/10.5281/zenodo.7057076)
Example output data from Fractal in the OME-Zarr format can be found here: [10.5281/zenodo.7081622](https://doi.org/10.5281/zenodo.7081622)
Example workflows can be found in the `examples` folder, together with additional instructions for how to set up the server & client, download the test data and run workflows through Fractal.

Fractal is currently in an early alpha version. We have the core processing functionality working for Yokogawa CV7000 image data and a workflow for processing OME-Zarr images up to feature measurements. But we're still adding core functionality and will introduce breaking changes. You can follow along our planned milestones on the [architecture side](https://github.com/fractal-analytics-platform/fractal/milestones) & the [tasks side](https://github.com/fractal-analytics-platform/fractal-tasks-core). Open an issue to get in touch, raise bugs or ask questions.

OME-Zarr files can be interactively visualizated in napari. Here is an example using the newly-proposed async loading in [NAP4](https://github.com/napari/napari/pull/4905) and the [napari-ome-zarr plugin](https://github.com/ome/napari-ome-zarr):

![napari_plate_overview](https://user-images.githubusercontent.com/18033446/190983839-afb9743f-530c-4b00-bde7-23ad62404ee8.gif)

### Contributors
Fractal was conceived in the Liberali Lab at the Friedrich Miescher Institute for Biomedical Research and in the Pelkmans Lab at the University of Zurich (both in Switzerland). The project lead is with [@gusqgm](https://github.com/gusqgm) & [@jluethi](https://github.com/jluethi).
The core development is done under contract by [@mfranzon](https://github.com/mfranzon), [@tcompa](https://github.com/tcompa) & [jacopo-exact](https://github.com/jacopo-exact) from eXact lab S.r.l. <exact-lab.it>.

## Installation

Simply

``` pip install fractal-client ```

Subsequently, you may invoke it as `fractal`. Note that you must provide
the following environment variables:

* `FRACTAL_SERVER`: fully qualified URL to the Fractal server installation
* `FRACTAL_USER`, `FRACTAL_PASSWORD`: email and password used to log-in to the
   Fractal server

By default, `fractal` caches some information in `~/.fractal`. This destination
can be customized by setting `SESSION_CACHE_PATH`.

For ease of use, you may define an environment file `.fractal.env` in the
folder from which `fractal` is invoked.


## Development

Development takes place on Github. You are welcome to submit an issue and open
pull requests.

### Developmente installation

Fractal is developed and maintained using [poetry](https://python-poetry.org/).

After cloning the repo, use

```
poetry install --with dev
```

to set up the development environment and all the dependencies and
dev-dependencies. You may run the test suite with

```
poetry run pytest
```

### Releasing

Before release checklist:

- [ ] The `main` branch is checked out
- [ ] You reviewed dependencies and dev dependencies and the lock file is up to
      date with `pyproject.toml`.
- [ ] The current `HEAD` of the main branch passes all the tests
- [ ] Use
```
poetry run bumpver update --dry --[patch|minor] --tag-commit --commit
```
to test updating the version bump
- [ ] If the previous step looks good, use
```
poetry run bumpver update --[patch|minor] --tag-commit --commit
```
to actually bump the version and commit the changes locally.
- [ ] Test the build with
```
poetry build
```
- [ ] If the previous step was successful, push the version bump and tags:
```
git push && git push --tags
```
- [ ] Finally, publish the updated package to pypi with
```
poetry publish --dry-run
```
removing the `--dry-run` when you made sure that everything looks good.
