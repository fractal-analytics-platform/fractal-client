"""
Copyright 2022 (C) Friedrich Miescher Institute for Biomedical Research and
University of Zurich

Original authors:
Jacopo Nespolo <jacopo.nespolo@exact-lab.it>
Marco Franzon <marco.franzon@exact-lab.it>
Tommaso Comparin <tommaso.comparin@exact-lab.it>

This file is part of Fractal and was originally developed by eXact lab S.r.l.
<exact-lab.it> under contract with Liberali Lab from the Friedrich Miescher
Institute for Biomedical Research and Pelkmans Lab from the University of
Zurich.
"""
import argparse as ap


parser_main = ap.ArgumentParser(description="Fractal Analytics Framework")

parser_main.add_argument(
    "-u",
    "--user",
    help="User email address for login (overrides configuration file)",
)
parser_main.add_argument(
    "-p", "--password", help="User password (overrides cnofiguration file)"
)
parser_main.add_argument(
    "-v",
    action="count",
    default=0,
    help="Use one or more time to increase verbosity level",
)
parser_main.add_argument(
    "--batch",
    default=False,
    action="store_true",
    help=(
        "Return output suitable for scripting, e.g., "
        "only the id of items created instead of the full object."
    ),
)

subparsers_main = parser_main.add_subparsers(title="Commands:", dest="cmd")

# REGISTER GROUP
project_parser = subparsers_main.add_parser(
    "register", help="Register with the Fractal server"
)

# PROJECT GROUP
project_parser = subparsers_main.add_parser("project", help="project commands")
project_subparsers = project_parser.add_subparsers(
    title="Valid subcommands:", dest="subcmd"
)

# PROJECT NEW
project_new_parser = project_subparsers.add_parser(
    "new", help="Create new project"
)
project_new_parser.add_argument("name", help="Name of new project")
project_new_parser.add_argument(
    "path",
    help=(
        "Project directory of new project. "
        "New datasets will be written here."
    ),
)
project_new_parser.add_argument(
    "-d",
    "--dataset",
    help=(
        "Name of new dataset to create. "
        "By default, the dataset `default` is created."
    ),
)

project_delete_parser = project_subparsers.add_parser(
    "delete", help="Create new project"
)

# DATASET GROUP
project_parser = subparsers_main.add_parser("dataset", help="dataset commands")

# TASK GROUP
project_parser = subparsers_main.add_parser("task", help="task commands")

# VERSION GROUP
project_parser = subparsers_main.add_parser(
    "version", help="Print verison and exit"
)


if __name__ == "__main__":
    args = parser_main.parse_args()

    from devtools import debug

    debug(args)
