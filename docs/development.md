# Contribute to Fractal Client development

The development of Fractal Client takes place on the [fractal-client Github
repository](https://github.com/fractal-analytics-platform/fractal-client).  To
ask questions or to inform us of a bug or unexpected behavior, please feel free
to [open an issue](https://github.com/fractal-analytics-platform/fractal-client/issues/new).


## Set up the development environment

### Clone repository

First, you should clone the repository
```
git clone https://github.com/fractal-analytics-platform/fractal-client.git
cd fractal-client
```

### Install package

We use [uv](https://docs.astral.sh/uv/) to manage the development environment and the dependencies - see https://docs.astral.sh/uv/getting-started/installation/ for methods to install it.
Running
```console
$ uv venv
$ uv sync --frozen [--group dev] [--group docs]
```
will create a new virtual environment in `./.venv` and install the main dependencies (and optionally the dev/docs groups).


## Build and release

Preliminary check-list:

* The `main` branch is checked out.
* You reviewed dependencies, and the lock file is up to date with `pyproject.toml`.
* The current HEAD of the `main` branch passes all the tests (`uv run pytest`).
* You updated the `CHANGELOG.md` file.
* You updated `docs/versions.md` with the constraints for the new version.

Actual **release instructions**:

1. Use one of the following
```
uv run --frozen bumpver update --tag-num --tag-commit --commit --dry
uv run --frozen bumpver update --patch --tag-commit --commit --dry
uv run --frozen bumpver update --minor --tag-commit --commit --dry
uv run --frozen bumpver update --set-version X.Y.Z --tag-commit --commit --dry
```
to test updating the version bump.

2. If the previous step looks good, remove `--dry` and re-run to actually bump the
version. This will trigger a dedicated GitHub action to build the new package
and publish it to PyPI.


## Run tests

Unit and integration testing of Fractal Server uses the
[pytest](https://docs.pytest.org/en/7.1.x/) testing framework.

Before running the test suite, a Fractal Server instance and a PostgreSQL
database must be available. For convenience, both services can be started via
Docker Compose using the configuration provided at `tests/fractal-server/docker-compose.yml`:
```sh
docker compose -f tests/fractal-server/docker-compose.yml up
```

The Fractal Server service can be customized via environment variables:

- `FRACTAL_SERVER_PYTHON`: Python version used to run Fractal Server (defaults to `3.12`);
- `FRACTAL_SERVER_BRANCH`: Git branch of the Fractal Server repository to use (defaults to `main`).

For example:
```sh
export FRACTAL_SERVER_PYTHON=3.14
export FRACTAL_SERVER_BRANCH=123-bug-fix
docker compose -f tests/fractal-server/docker-compose.yml up
```

Then, you may run the test suite by invoking
```sh
uv run pytest
```
from the main directory of the `fractal-client` repository. It is sometimes
useful to specify additional arguments, e.g.
```sh
uv run pytest -s -vvv --log-cli-level info --full-trace
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
uv run mkdocs serve --config-file mkdocs.yml  # serves the docs at http://127.0.0.1:8000
uv run mkdocs build --config-file mkdocs.yml  # creates a build in the `site` folder
```
