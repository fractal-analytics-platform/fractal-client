# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/fractal-analytics-platform/fractal/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                                       |    Stmts |     Miss |   Branch |   BrPart |   Cover |   Missing |
|------------------------------------------- | -------: | -------: | -------: | -------: | ------: | --------: |
| fractal/\_\_init\_\_.py                    |        1 |        0 |        0 |        0 |    100% |           |
| fractal/authclient.py                      |       70 |        5 |       12 |        3 |     90% |27->exit, 28, 43-44, 75-76 |
| fractal/client.py                          |       56 |       10 |       18 |        4 |     81% |60-67, 76, 86-87, 112, 118-120, 124 |
| fractal/cmd/\_\_init\_\_.py                |      224 |        6 |       82 |        6 |     96% |93, 141, 203, 281, 308, 368 |
| fractal/cmd/\_aux\_task\_caching.py        |       84 |        0 |       35 |        0 |    100% |           |
| fractal/cmd/\_dataset.py                   |       97 |        2 |       28 |        2 |     97% |   90, 142 |
| fractal/cmd/\_job.py                       |       73 |       11 |       22 |        2 |     86% |106-114, 120-131, 163 |
| fractal/cmd/\_project.py                   |       67 |        4 |       22 |        3 |     92% |28, 38-42, 66 |
| fractal/cmd/\_task.py                      |      119 |        4 |       64 |        4 |     96% |69-70, 123->125, 124, 202 |
| fractal/cmd/\_user.py                      |       73 |        6 |       30 |        2 |     90% |41-46, 121 |
| fractal/cmd/\_workflow.py                  |      133 |        7 |       54 |        7 |     91% |42, 155, 165->170, 240, 286->294, 295-298, 327->335 |
| fractal/common/\_\_init\_\_.py             |        0 |        0 |        0 |        0 |    100% |           |
| fractal/common/schemas/\_\_init\_\_.py     |        9 |        0 |        0 |        0 |    100% |           |
| fractal/common/schemas/\_validators.py     |       31 |        1 |       16 |        1 |     96% |        54 |
| fractal/common/schemas/applyworkflow.py    |       49 |        0 |       12 |        0 |    100% |           |
| fractal/common/schemas/manifest.py         |       41 |        0 |       12 |        0 |    100% |           |
| fractal/common/schemas/project.py          |       50 |        0 |        0 |        0 |    100% |           |
| fractal/common/schemas/state.py            |       14 |        0 |        0 |        0 |    100% |           |
| fractal/common/schemas/task.py             |       62 |        0 |        0 |        0 |    100% |           |
| fractal/common/schemas/task\_collection.py |       39 |        4 |        8 |        1 |     89% |74, 104-106 |
| fractal/common/schemas/user.py             |       24 |        0 |        0 |        0 |    100% |           |
| fractal/common/schemas/workflow.py         |       56 |        0 |       11 |        0 |    100% |           |
| fractal/config.py                          |       22 |        4 |        2 |        0 |     75% |     25-28 |
| fractal/interface.py                       |       37 |        1 |        4 |        0 |     98% |        18 |
| fractal/parser.py                          |      193 |        0 |        0 |        0 |    100% |           |
| fractal/response.py                        |       22 |        0 |        6 |        0 |    100% |           |
|                                  **TOTAL** | **1646** |   **65** |  **438** |   **35** | **95%** |           |


## Setup coverage badge

Below are examples of the badges you can use in your main branch `README` file.

### Direct image

[![Coverage badge](https://raw.githubusercontent.com/fractal-analytics-platform/fractal/python-coverage-comment-action-data/badge.svg)](https://htmlpreview.github.io/?https://github.com/fractal-analytics-platform/fractal/blob/python-coverage-comment-action-data/htmlcov/index.html)

This is the one to use if your repository is private or if you don't want to customize anything.

### [Shields.io](https://shields.io) Json Endpoint

[![Coverage badge](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/fractal-analytics-platform/fractal/python-coverage-comment-action-data/endpoint.json)](https://htmlpreview.github.io/?https://github.com/fractal-analytics-platform/fractal/blob/python-coverage-comment-action-data/htmlcov/index.html)

Using this one will allow you to [customize](https://shields.io/endpoint) the look of your badge.
It won't work with private repositories. It won't be refreshed more than once per five minutes.

### [Shields.io](https://shields.io) Dynamic Badge

[![Coverage badge](https://img.shields.io/badge/dynamic/json?color=brightgreen&label=coverage&query=%24.message&url=https%3A%2F%2Fraw.githubusercontent.com%2Ffractal-analytics-platform%2Ffractal%2Fpython-coverage-comment-action-data%2Fendpoint.json)](https://htmlpreview.github.io/?https://github.com/fractal-analytics-platform/fractal/blob/python-coverage-comment-action-data/htmlcov/index.html)

This one will always be the same color. It won't work for private repos. I'm not even sure why we included it.

## What is that?

This branch is part of the
[python-coverage-comment-action](https://github.com/marketplace/actions/python-coverage-comment)
GitHub Action. All the files in this branch are automatically generated and may be
overwritten at any moment.