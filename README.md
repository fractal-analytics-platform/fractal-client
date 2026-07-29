# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/fractal-analytics-platform/fractal-client/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                                            |    Stmts |     Miss |   Branch |   BrPart |   Cover |   Missing |
|------------------------------------------------ | -------: | -------: | -------: | -------: | ------: | --------: |
| src/fractal\_client/\_\_init\_\_.py             |        1 |        0 |        0 |        0 |    100% |           |
| src/fractal\_client/auth/\_\_init\_\_.py        |        0 |        0 |        0 |        0 |    100% |           |
| src/fractal\_client/auth/\_args.py              |       24 |        0 |        4 |        0 |    100% |           |
| src/fractal\_client/auth/\_client.py            |       64 |        0 |        6 |        0 |    100% |           |
| src/fractal\_client/auth/\_cmd\_handler.py      |       11 |        0 |        4 |        0 |    100% |           |
| src/fractal\_client/auth/\_token\_utils.py      |       40 |        0 |        8 |        0 |    100% |           |
| src/fractal\_client/cmd/\_\_init\_\_.py         |      292 |        0 |      100 |        0 |    100% |           |
| src/fractal\_client/cmd/\_aux\_task\_caching.py |       83 |        0 |       20 |        0 |    100% |           |
| src/fractal\_client/cmd/\_dataset.py            |       25 |        0 |        6 |        0 |    100% |           |
| src/fractal\_client/cmd/\_group.py              |       28 |        0 |        4 |        0 |    100% |           |
| src/fractal\_client/cmd/\_job.py                |       72 |        1 |       22 |        0 |     99% |       104 |
| src/fractal\_client/cmd/\_profile.py            |       13 |        0 |        2 |        0 |    100% |           |
| src/fractal\_client/cmd/\_project.py            |       29 |        0 |        4 |        0 |    100% |           |
| src/fractal\_client/cmd/\_resource.py           |       13 |        0 |        2 |        0 |    100% |           |
| src/fractal\_client/cmd/\_task.py               |       72 |        0 |       30 |        0 |    100% |           |
| src/fractal\_client/cmd/\_task\_collection.py   |       62 |        0 |       26 |        0 |    100% |           |
| src/fractal\_client/cmd/\_template.py           |       44 |        0 |       14 |        0 |    100% |           |
| src/fractal\_client/cmd/\_user.py               |       70 |        0 |       30 |        1 |     99% |  99-\>101 |
| src/fractal\_client/cmd/\_workflow.py           |      143 |        0 |       48 |        2 |     99% |215-\>214, 225-\>227 |
| src/fractal\_client/config.py                   |       17 |        0 |        0 |        0 |    100% |           |
| src/fractal\_client/interface.py                |       10 |        0 |        2 |        0 |    100% |           |
| src/fractal\_client/main.py                     |       43 |        0 |        6 |        0 |    100% |           |
| src/fractal\_client/parser/\_\_init\_\_.py      |       25 |        0 |        0 |        0 |    100% |           |
| src/fractal\_client/parser/\_dataset.py         |       14 |        0 |        0 |        0 |    100% |           |
| src/fractal\_client/parser/\_job.py             |       25 |        0 |        0 |        0 |    100% |           |
| src/fractal\_client/parser/\_main.py            |       11 |        0 |        0 |        0 |    100% |           |
| src/fractal\_client/parser/\_profile.py         |        6 |        0 |        0 |        0 |    100% |           |
| src/fractal\_client/parser/\_project.py         |       18 |        0 |        0 |        0 |    100% |           |
| src/fractal\_client/parser/\_resource.py        |        5 |        0 |        0 |        0 |    100% |           |
| src/fractal\_client/parser/\_task.py            |       50 |        0 |        0 |        0 |    100% |           |
| src/fractal\_client/parser/\_template.py        |       17 |        0 |        0 |        0 |    100% |           |
| src/fractal\_client/parser/\_user.py            |       29 |        0 |        0 |        0 |    100% |           |
| src/fractal\_client/parser/\_usergroup.py       |       15 |        0 |        0 |        0 |    100% |           |
| src/fractal\_client/parser/\_version.py         |        2 |        0 |        0 |        0 |    100% |           |
| src/fractal\_client/parser/\_workflow.py        |       56 |        0 |        0 |        0 |    100% |           |
| src/fractal\_client/response.py                 |       40 |        0 |       12 |        0 |    100% |           |
| **TOTAL**                                       | **1469** |    **1** |  **350** |    **3** | **99%** |           |


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