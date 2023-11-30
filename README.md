# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/fractal-analytics-platform/fractal-client/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                                        |    Stmts |     Miss |   Branch |   BrPart |   Cover |   Missing |
|-------------------------------------------- | -------: | -------: | -------: | -------: | ------: | --------: |
| fractal\_client/\_\_init\_\_.py             |        1 |        0 |        0 |        0 |    100% |           |
| fractal\_client/authclient.py               |       70 |        4 |       12 |        1 |     94% |43-44, 75-76 |
| fractal\_client/client.py                   |       56 |       10 |       18 |        4 |     81% |60-67, 76, 86-87, 112, 118-120, 124 |
| fractal\_client/cmd/\_\_init\_\_.py         |      219 |        6 |       80 |        6 |     96% |98, 146, 208, 286, 313, 368 |
| fractal\_client/cmd/\_aux\_task\_caching.py |       84 |        0 |       35 |        0 |    100% |           |
| fractal\_client/cmd/\_dataset.py            |       83 |        1 |       26 |        1 |     98% |       122 |
| fractal\_client/cmd/\_job.py                |       70 |       11 |       20 |        2 |     86% |99-107, 113-124, 156 |
| fractal\_client/cmd/\_project.py            |       51 |        1 |       14 |        1 |     97% |        45 |
| fractal\_client/cmd/\_task.py               |      109 |        1 |       64 |        1 |     99% |       185 |
| fractal\_client/cmd/\_user.py               |       66 |        6 |       32 |        2 |     90% |37-42, 105 |
| fractal\_client/cmd/\_workflow.py           |      121 |        1 |       56 |        4 |     97% |139, 149->154, 263->271, 302->310 |
| fractal\_client/config.py                   |       17 |        0 |        2 |        0 |    100% |           |
| fractal\_client/interface.py                |       37 |        1 |        4 |        0 |     98% |        18 |
| fractal\_client/parser.py                   |      191 |        0 |        0 |        0 |    100% |           |
| fractal\_client/response.py                 |       29 |        0 |        6 |        0 |    100% |           |
|                                   **TOTAL** | **1204** |   **42** |  **369** |   **22** | **96%** |           |


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