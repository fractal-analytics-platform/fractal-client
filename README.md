# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/fractal-analytics-platform/fractal-client/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                                        |    Stmts |     Miss |   Branch |   BrPart |   Cover |   Missing |
|-------------------------------------------- | -------: | -------: | -------: | -------: | ------: | --------: |
| fractal\_client/\_\_init\_\_.py             |        1 |        0 |        0 |        0 |    100% |           |
| fractal\_client/authclient.py               |       78 |        3 |        6 |        0 |     96% | 38, 82-83 |
| fractal\_client/client.py                   |       57 |        4 |       12 |        1 |     93% |131-133, 137 |
| fractal\_client/cmd/\_\_init\_\_.py         |      243 |        3 |       84 |        1 |     99% |   289-291 |
| fractal\_client/cmd/\_aux\_task\_caching.py |       83 |        1 |       20 |        1 |     98% |       231 |
| fractal\_client/cmd/\_dataset.py            |       27 |        0 |        6 |        1 |     97% |    45->48 |
| fractal\_client/cmd/\_group.py              |       36 |        0 |        6 |        0 |    100% |           |
| fractal\_client/cmd/\_job.py                |       73 |        8 |       22 |        1 |     91% |82-93, 117-121 |
| fractal\_client/cmd/\_project.py            |       30 |        0 |        4 |        0 |    100% |           |
| fractal\_client/cmd/\_task.py               |       65 |        0 |       30 |        0 |    100% |           |
| fractal\_client/cmd/\_task\_collection.py   |       60 |        0 |       26 |        0 |    100% |           |
| fractal\_client/cmd/\_user.py               |      120 |        0 |       50 |        1 |     99% |    75->84 |
| fractal\_client/cmd/\_workflow.py           |      120 |        1 |       36 |        1 |     99% |        96 |
| fractal\_client/config.py                   |       14 |        0 |        0 |        0 |    100% |           |
| fractal\_client/interface.py                |       10 |        0 |        2 |        0 |    100% |           |
| fractal\_client/parser.py                   |      203 |        0 |        0 |        0 |    100% |           |
| fractal\_client/response.py                 |       40 |        1 |       12 |        1 |     96% |        42 |
|                                   **TOTAL** | **1260** |   **21** |  **316** |    **8** | **98%** |           |


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