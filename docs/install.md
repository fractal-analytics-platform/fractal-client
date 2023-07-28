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
command line (see its documentation [here](/reference/fractal/)).  Note that you must provide
the following environment variables:

* `FRACTAL_SERVER`: fully qualified URL to the Fractal server installation
* `FRACTAL_USER`, `FRACTAL_PASSWORD`: email and password used to log-in to the
   Fractal server.

By default, `fractal` caches some information in `~/.cache/fractal`. This destination
can be customized by setting `FRACTAL_CACHE_PATH`.

For ease of use, you may define an environment file `.fractal.env` in the
folder from which `fractal` is invoked, with the relevant environment
variables.  An example of such file is
```
FRACTAL_USER=user@something.com
FRACTAL_PASSWORD=myuser
FRACTAL_SERVER=http://localhost:8010
```
