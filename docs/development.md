# Contribute to Fractal Client development

The development of Fractal Client takes place on the [fractal-client Github
repository](https://github.com/fractal-analytics-platform/fractal-client).  To
ask questions or to inform us of a bug or unexpected behavior, please feel free
to [open an
issue](https://github.com/fractal-analytics-platform/fractal-client/issues/new).


## Set up the development environment

### Clone repository

First, you should clone the repository
```
git clone https://github.com/fractal-analytics-platform/fractal-client.git
cd fractal
```

### Install package

We use [poetry](https://python-poetry.org/docs) v1.5 to manage the development environment and the dependencies. A simple way to install it is `pipx install poetry==1.7.1`, or you can look at the installation section [here](https://python-poetry.org/docs#installation).
Running

```console
poetry install [--with dev] [--with docs]
```
will take care of installing all the dependencies in a separate environment, optionally installing also the dependencies for developement and to build the documentation.


## Build and release

We also use `poetry` to build the package and publish it to PyPI.

Preliminary check-list:

* The `main` branch is checked out.
* You reviewed dependencies, and the lock file is up to date with ``pyproject.toml``.
* The current HEAD of the `main` branch passes all the tests (note: make sure
  that you are using `poetry run pytest`, and not simply `pytest`).
* You updated the `CHANGELOG.md` file.

Actual **release instructions**:

1. Use:
```
poetry run bumpver update --[tag-num|patch|minor] --tag-commit --commit --dry
```
to test updating the version bump.

2. If the previous step looks good, remove `--dry` and re-run to actually bump the
version. This will trigger a dedicated GitHub action to build the new package
and publish it to PyPI.


## Run tests

Unit and integration testing of Fractal Server uses the
[pytest](https://docs.pytest.org/en/7.1.x/) testing framework.

If you installed the development dependencies, you may run
the test suite by invoking
```
poetry run pytest
```
from the main directory of the `fractal-client` repository. It is sometimes
useful to specify additional arguments, e.g.
```
poetry run pytest -s -vvv --log-cli-level info --full-trace
```

Tests are also run as part of [GitHub Actions Continuous
Integration](https://github.com/fractal-analytics-platform/fractal-client/actions/workflows/ci.yml)
for the `fractal-client` repository.


## Documentation

The documentations is built with mkdocs, and we bundle a module from
[sphinx-argparse plugin](https://sphinx-argparse.readthedocs.io), customized to
our needs.

To build the documentation locally, setup a development python environment (e.g. with `poetry install --with docs`) and then run one of these commands:
```
poetry run mkdocs serve --config-file mkdocs.yml  # serves the docs at http://127.0.0.1:8000
poetry run mkdocs build --config-file mkdocs.yml  # creates a build in the `site` folder
```
