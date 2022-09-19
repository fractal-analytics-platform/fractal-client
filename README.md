# Fractal Client

The client component of Fractal Analytics Framework.

[![PyPI
version](https://img.shields.io/pypi/v/fractal-client?color=gree)](https://pypi.org/project/fractal-client/)

Fractal is a framework to process high content screening data at scale and
prepares it for interactive visualization.

Fractal provides distributed workflows that convert TBs of image data into
OME-Zarr files. The platform then processes the 3D image data by applying tasks
like illumination correction and maximum intensity projection. The pyramidal
OME-Zarr files enable interactive visualization in the napari viewer. We are
building towards integrating object segmentation (nuclei, cells, organoids) and
feature measurements into Fractal.

Fractal is currently an **early alpha build**. We currently support only
Yokogawa CV7000 image data as an input. Also, we're in the process of
refactoring the workflow management into a client-server architecture. Thus,
proceed with it at your own risk, there will still be significant breaking
changes. While we don't have any releases or stable versions and thus don't
provide user support, we welcome questions and discussions. Open an issue to
get in touch.


![Fractal_multiwell_plate](https://user-images.githubusercontent.com/18033446/177169496-09784413-6ba9-4041-80e2-c70a70e0a5d9.gif)

Shortened movie of browsing an OME-Zarr file generated with Fractal in napari,
using the [napari-ome-zarr plugin](https://github.com/ome/napari-ome-zarr).
Actual loading times vary and can be a bit slower than in this GIF.

### Contributors

Fractal was conceived in the Liberali Lab at the Friedrich Miescher Institute
for Biomedical Research and in the Pelkmans Lab at the University of Zurich
(both in Switzerland). The project lead is with
[@gusqgm](https://github.com/gusqgm) & [@jluethi](https://github.com/jluethi).
The project was originally led by [@dvischi](https://github.com/dvischi). The
core development is done under contract by
[@mfranzon](https://github.com/mfranzon), [@tcompa](https://github.com/tcompa)
& [jacopo-exact](https://github.com/jacopo-exact) from eXact lab S.r.l.
<exact-lab.it>.


## Installation

Simply

``` pip install fractal-client ```

Subsequently, you may simply invoke it as `fractal`. Note that you must provide
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
