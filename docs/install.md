# Installation and usage

## Installation

Fractal Client is hosted on [the PyPI
index](https://pypi.org/project/fractal-client), and it can be installed with
`pip` via
```
pip install fractal-client
```

## Usage

You may invoke the Fractal Client via the custom command `fractal`, from the
command line (see its documentation [here](/reference/fractal/)).

You must set the `FRACTAL_SERVER` variable, which is a fully qualified URL to
the Fractal server installation (e.g. http://localhost:8000). This can be an
environment variable or it can be stored in a an environment file
`.fractal.env` as in
```
FRACTAL_SERVER=http://localhost:8010
```

### Credentials

Most `fractal` commands are restricted to authenticated users, and user
credentials can be specified in multiple ways:
* Set `FRACTAL_USER` and `FRACTAL_PASSWORD` variables as environment variables;
* Add `FRACTAL_USER` and `FRACTAL_PASSWORD` variables in `.fractal.env`;
* Explicitly provide `--user` and `--password` arguments for `fractal` commands, see [here](/reference/fractal/).

### Cache

By default, `fractal` caches some information (namely a valid token for the
current session on `fractal-server` and a list of tasks) in `~/.cache/fractal`.
This destination can be customized by setting the `FRACTAL_CACHE_PATH`
variables.

### Full example

Here is an example of a valid `.fractal.env` file:
```
FRACTAL_USER=user@something.com
FRACTAL_PASSWORD=myuser
FRACTAL_SERVER=http://localhost:8010
FRACTAL_CACHE_PATH=/some/path/fractal-cache
```
