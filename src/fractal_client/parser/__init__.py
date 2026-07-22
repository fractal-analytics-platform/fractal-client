"""
Copyright 2022 (C) Friedrich Miescher Institute for Biomedical Research and
University of Zurich

Original authors:
Jacopo Nespolo <jacopo.nespolo@exact-lab.it>
Marco Franzon <marco.franzon@exact-lab.it>
Tommaso Comparin <tommaso.comparin@exact-lab.it>
Yuri Chiucconi  <yuri.chiucconi@exact-lab.it>

This file is part of Fractal and was originally developed by eXact lab S.r.l.
<exact-lab.it> under contract with Liberali Lab from the Friedrich Miescher
Institute for Biomedical Research and Pelkmans Lab from the University of
Zurich.
"""

from ._dataset import add_dataset_parser
from ._job import add_job_parser
from ._main import get_main_parser
from ._profile import add_profile_parser
from ._project import add_project_parser
from ._resource import add_resource_parser
from ._task import add_task_parser
from ._template import add_template_parser
from ._user import add_user_parser
from ._version import add_version_parser
from ._workflow import add_workflow_parser

parser_main = get_main_parser()

subparsers_main = parser_main.add_subparsers(title="Commands", dest="cmd")
add_project_parser(subparsers_main)
add_dataset_parser(subparsers_main)
add_task_parser(subparsers_main)
add_workflow_parser(subparsers_main)
add_job_parser(subparsers_main)
add_version_parser(subparsers_main)
add_user_parser(subparsers_main)
add_resource_parser(subparsers_main)
add_profile_parser(subparsers_main)
add_template_parser(subparsers_main)
