# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/fractal-analytics-platform/fractal-client/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                                        |    Stmts |     Miss |   Branch |   BrPart |   Cover |   Missing |
|-------------------------------------------- | -------: | -------: | -------: | -------: | ------: | --------: |
| fractal\_client/\_\_init\_\_.py             |        1 |        0 |        0 |        0 |    100% |           |
| fractal\_client/authclient.py               |       78 |        2 |       14 |        1 |     97% |     82-83 |
| fractal\_client/client.py                   |       58 |        4 |       20 |        1 |     94% |132-134, 138 |
| fractal\_client/cmd/\_\_init\_\_.py         |      203 |        0 |       72 |        0 |    100% |           |
| fractal\_client/cmd/\_aux\_task\_caching.py |       84 |        0 |       27 |        0 |    100% |           |
| fractal\_client/cmd/\_dataset.py            |       35 |        0 |       12 |        0 |    100% |           |
| fractal\_client/cmd/\_job.py                |       67 |        6 |       24 |        1 |     92% |82-93, 121 |
| fractal\_client/cmd/\_project.py            |       31 |        0 |        4 |        0 |    100% |           |
| fractal\_client/cmd/\_task.py               |      125 |        0 |       76 |        0 |    100% |           |
| fractal\_client/cmd/\_user.py               |       65 |        0 |       34 |        1 |     99% |    38->47 |
| fractal\_client/cmd/\_workflow.py           |      127 |        0 |       62 |        2 |     99% |242->250, 281->289 |
| fractal\_client/config.py                   |       15 |        0 |        2 |        0 |    100% |           |
| fractal\_client/interface.py                |       10 |        0 |        2 |        0 |    100% |           |
| fractal\_client/parser.py                   |      183 |        0 |        0 |        0 |    100% |           |
| fractal\_client/response.py                 |       29 |        0 |        6 |        0 |    100% |           |
|                                   **TOTAL** | **1111** |   **12** |  **355** |    **6** | **99%** |           |


## Setup coverage badge

Below are examples of the badges you can use in your main branch `README` file.

### Direct image

[![Coverage badge](https://raw.githubusercontent.com/fractal-analytics-platform/fractal-client/python-coverage-comment-action-data/badge.svg)](https://htmlpreview.github.io/?https://github.com/fractal-analytics-platform/fractal-client/blob/python-coverage-comment-action-data/htmlcov/index.html)

This is the one to use if your repository is private or if you don't want to customize anything.

### [Shields.io](https://shields.io) Json Endpoint

[![Coverage badge](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/fractal-analytics-platform/fractal-client/python-coverage-comment-action-data/endpoint.json)](https://htmlpreview.github.io/?https://github.com/fractal-analytics-platform/fractal-client/blob/python-coverage-comment-action-data/htmlcov/index.html)

Using this one will allow you to [customize](https://shields.io/endpoint) the look of your badge.
It won't work with private repositories. It won't be refreshed more than once per five minutes.

### [Shields.io](https://shields.io) Dynamic Badge

[![Coverage badge](https://img.shields.io/badge/dynamic/json?color=brightgreen&label=coverage&query=%24.message&url=https%3A%2F%2Fraw.githubusercontent.com%2Ffractal-analytics-platform%2Ffractal-client%2Fpython-coverage-comment-action-data%2Fendpoint.json)](https://htmlpreview.github.io/?https://github.com/fractal-analytics-platform/fractal-client/blob/python-coverage-comment-action-data/htmlcov/index.html)

This one will always be the same color. It won't work for private repos. I'm not even sure why we included it.

## What is that?

This branch is part of the
[python-coverage-comment-action](https://github.com/marketplace/actions/python-coverage-comment)
GitHub Action. All the files in this branch are automatically generated and may be
overwritten at any moment.