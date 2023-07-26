# Contribute to Fractal Client development

The development of Fractal Client takes place on the [fractal Github
repository](https://github.com/fractal-analytics-platform/fractal).  To
ask questions or to inform us of a bug or unexpected behavior, please feel free
to [open an
issue](https://github.com/fractal-analytics-platform/fractal/issues/new).

To contribute code, please fork the repository and submit a pull request.

## Set up the development environment

### Install poetry

Fractal uses [poetry](https://python-poetry.org/docs) to manage the development
environment and dependencies, and to streamline the build and release
operations. Version 1.3 is recommended, although 1.2.2 should also work.

A simple way to install it is the command `pipx install poetry==1.3`; other
options are described
[here](https://python-poetry.org/docs#installing-with-the-official-installer).


### Clone repositories

You can clone the `fractal` repository via
```
git clone https://github.com/fractal-analytics-platform/fractal.git
```
and then (from the `fractal` folder)
```
git submodule update --init
```
The second command is needed, since `fractal` includes
[`fractal-common`](https://github.com/fractal-analytics-platform/fractal-common)
as a git submodule.

### Install package

Running
```
poetry install
```
will initialise a Python virtual environment and install Fractal Server and all
its dependencies, including development dependencies.
Note that to run commands from within this environment you should prepend them
with `poetry run`, as in `poetry run fractalctl start`.

It may sometimes be useful to use a different Python interpreter from the one
installed in your system. To this end we suggest using
[pyenv](https://github.com/pyenv/pyenv). In the project folder, invoking
```
pyenv local <version>
poetry env use <version>
```
will install Fractal in a development environment using an interpreter pinned
at the version provided instead of the system interpreter.

## Build and release

Preliminary check-list

* The `main` branch is checked out.
* You reviewed dependencies, and the lock file is up to date with ``pyproject.toml``.
* The current HEAD of the `main` branch passes all the tests (note: make sure
  that you are using the poetry-installed local package).
* You updated the `CHANGELOG.md` file.

Actual **release instructions**:

1. Use:
```
poetry run bumpver update --[tag-num|patch|minor] --tag-commit --commit --dry
```
to test updating the version bump.

2. If the previous step looks good, remove `--dry` and re-run to actually bump the
version and commit the changes locally.

3. Test the build with:
```
poetry build
```
4. If the previous step was successful, push the version bump and tags:
```
git push && git push --tags
```
5. Finally, publish the updated package to PyPI with:
```
poetry publish --dry-run
```
replacing ``--dry-run`` with ``--username YOUR_USERNAME --password
YOUR_PASSWORD`` when you made sure that everything looks good.


## Run tests

Unit and integration testing of Fractal Server uses the
[pytest](https://docs.pytest.org/en/7.1.x/) testing framework.

If you installed the development dependencies, you may run
the test suite by invoking
```
poetry run pytest
```
from the main directory of the `fractal` repository. It is sometimes
useful to specify additional arguments, e.g.
```
poetry run pytest -s -vvv --log-cli-level info --full-trace
```

Tests are also run as part of [GitHub Actions Continuous
Integration](https://github.com/fractal-analytics-platform/fractal/actions/workflows/ci.yml)
for the `fractal` repository.


## Documentation

The documentations is built with [sphinx](https://www.sphinx-doc.org), the
[Read the Docs theme](https://sphinx-rtd-theme.readthedocs.io) and the
[sphinx-argparse plugin](https://sphinx-argparse.readthedocs.io).

To build the documentation locally

1. Setup a development python environment (e.g. with `poetry install --with docs`)
2. Run `poetry run sphinx-build docs/source docs/build`.
