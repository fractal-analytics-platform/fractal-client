# Fractal Client

[![PyPI version](https://img.shields.io/pypi/v/fractal-client?color=gree)](https://pypi.org/project/fractal-client/)
[![CI Status](https://github.com/fractal-analytics-platform/fractal/actions/workflows/ci.yml/badge.svg)](https://github.com/fractal-analytics-platform/fractal/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

Fractal is a framework to process high-content imaging data at scale and prepare it for interactive visualization.

![Fractal_Overview](https://fractal-analytics-platform.github.io/assets/fractal_overview.jpg)

Fractal provides distributed workflows that convert TBs of image data into OME-Zarr files. The platform then processes the 3D image data by applying tasks like illumination correction, maximum intensity projection, 3D segmentation using [cellpose](https://cellpose.readthedocs.io/en/latest/) and measurements using [napari workflows](https://github.com/haesleinhuepf/napari-workflows). The pyramidal OME-Zarr files enable interactive visualization in the napari viewer.

This is the repository that contains the **Fractal client**. Find more information about Fractal in general and the other repositories at the [Fractal home page](https://fractal-analytics-platform.github.io).

## Documentation

Documentation is not yet available.

For installation instructions, see below.

## Installation

Simply

```
pip install fractal-client
```

Subsequently, you may invoke it as `fractal`. Note that you must provide
the following environment variables:

* `FRACTAL_SERVER`: fully qualified URL to the Fractal server installation
* `FRACTAL_USER`, `FRACTAL_PASSWORD`: email and password used to log-in to the
   Fractal server

By default, `fractal` caches some information in `~/.cache/fractal`. This destination
can be customized by setting `FRACTAL_CACHE_PATH`.

For ease of use, you may define an environment file `.fractal.env` in the
folder from which `fractal` is invoked.


## Development

### Installation

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

# Contributors and license

Unless otherwise stated in each individual module, all Fractal components are released according to a BSD 3-Clause License, and Copyright is with Friedrich Miescher Institute for Biomedical Research and University of Zurich.

Fractal was conceived in the Liberali Lab at the Friedrich Miescher Institute for Biomedical Research and in the Pelkmans Lab at the University of Zurich (both in Switzerland). The project lead is with [@gusqgm](https://github.com/gusqgm) & [@jluethi](https://github.com/jluethi). The core development is done under contract by [@mfranzon](https://github.com/mfranzon), [@tcompa](https://github.com/tcompa) & [@jacopo-exact](https://github.com/jacopo-exact) from [eXact lab S.r.l.](exact-lab.it).
