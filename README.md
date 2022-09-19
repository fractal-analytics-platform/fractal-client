# Fractal Server

[![PyPI version](https://badge.fury.io/py/fractal-server.svg)](https://badge.fury.io/py/fractal-server)

Fractal is a framework to process high content imaging data at scale and prepare it for interactive visualization.

This is the server component of the fractal analytics platform. If you are interested in the client component, please refer to the [main
repository](https://github.com/fractal-analytics-platform/fractal). If you are interested in the fractal tasks, please refer to [the tasks repository](https://github.com/fractal-analytics-platform/fractal-tasks-core).

## Installation

You may
`pip install fractal-server`

This will install the project and its dependencies. If you wish to also install
Fractal core tasks, use
```
pip install fractal-server[core-tasks]
```

### Environment and database

You will need to define some environment variables in order to use
`fractal-server`. For your convenience, you may simply copy
`template.fractal_server.env` to `.fractal_server.env`.

`fractal-server` requires a database to run. Once you set up the environment
variables you need to initialise the database by invoking

```
alembic upgrade head
```

NOTE: as `fractal-server` is still in pre-alpha the database schema is not yet
committed to the repository. As such you'll be required to first create a
schema revision with `alembic revision --autogenerate -m "[your commit
message]"`.

## Contributing

To contribute to the development of `fractal-server` you may fork and clone the
[repository](https://github.com/fractal-analytics-platform/fractal-server).

We use [poetry](https://python-poetry.org/docs/) (v1.2) to manage the
development environment and the dependencies. Running

```
poetry install [--with dev]
```

will take care of installing all the dependencies in a separate environment,
optionally installing also the development dependencies.

It may be useful to use a different interpreter from the one installed in your
system. To this end we recommend using
[pyenv](https://github.com/pyenv/pyenv). In the project folder, invoking

```
pyenv local 3.8.13
poetry env use 3.8
poetry install
```

will install Fractal in a development environment using `python 3.8.13` instad
of the system-wide interpreter.

### Testing

We use [pytest](https://docs.pytest.org/en/7.1.x/) for unit and integration
testing of Fractal. If you installed the development dependencies, you may run
the test suite by invoking

```
poetry run pytest
```

## Contributors

Fractal was conceived in the Liberali Lab at the Friedrich Miescher Institute
for Biomedical Research and in the Pelkmans Lab at the University of Zurich
(both in Switzerland). The project lead is with
[@gusqgm](https://github.com/gusqgm) & [@jluethi](https://github.com/jluethi).
The core development is done under contract by
[@mfranzon](https://github.com/mfranzon), [@tcompa](https://github.com/tcompa)
& [@jacopo-exact](https://github.com/jacopo-exact) from [eXact lab
S.r.l.](https://exact-lab.it).

## License

Fractal is released according to a BSD 3-Clause License. See `LICENSE`.
