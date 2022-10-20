"""
Copyright 2022 (C) Friedrich Miescher Institute for Biomedical Research and
University of Zurich

Original authors:
Jacopo Nespolo <jacopo.nespolo@exact-lab.it>
Tommaso Comparin <tommaso.comparin@exact-lab.it>
Marco Franzon <marco.franzon@exact-lab.it>

This file is part of Fractal and was originally developed by eXact lab S.r.l.
<exact-lab.it> under contract with Liberali Lab from the Friedrich Miescher
Institute for Biomedical Research and Pelkmans Lab from the University of
Zurich.
"""
import json
import logging
from pathlib import Path
from sys import stdout
from typing import Any
from typing import Dict
from typing import Iterable
from typing import List
from typing import Optional

from pydantic import BaseModel


logger = logging.getLogger(__name__)


def dummy_fail(
    *,
    input_paths: Iterable[Path],
    output_path: Path,
    metadata: Optional[Dict[str, Any]] = None,
    # arguments of this task
    error_message: str,
) -> Dict[str, Any]:
    """
    Dummy task

    This task appends to a json file the parameters it was called with, such
    that it is easy to parse the file in a test settings.

    Incidentally, this task defines the reference interface of a task.

    Arguments
    ---------
    input_paths (iterable of Path) :
        The paths to fetch data from
    output_path (Path) :
        The output path, pointing either to a file or to a directory in which
        the task will write its output files.
    metadata (Dict or None) :
        Optional metadata about the input the task may need

    any further argument (Any) :
        the arguments specific to the issue, i.e., `errro_message` in the
        present example

    Retrun
    ------
    metadata_update (Dict[str, Any]) :
        a dictionary that will update the metadata
    """
    logger.info("ENTERING dummy_fail task")
    raise ValueError(error_message)

    return {}


if __name__ == "__main__":
    from argparse import ArgumentParser

    class TaskArguments(BaseModel):
        """
        Wrap task arguments to ease marshalling

        This way we can automatically cast the input from command line onto
        the correct type required by the task.
        """

        input_paths: List[Path]
        output_path: Path
        metadata: Optional[Dict[str, Any]] = None
        error_message: str

    parser = ArgumentParser()
    parser.add_argument("-j", "--json", help="Read parameters from json file")
    parser.add_argument(
        "-o",
        "--output",
        help=(
            "Output file to redirect serialised returned data "
            "(default stdout)"
        ),
    )

    args = parser.parse_args()

    if args.output and Path(args.output).exists():
        logger.error(f"Output file {args.output} already exists. Terminating")
        exit(1)

    pars = {}
    if args.json:
        with open(args.json, "r") as f:
            pars = json.load(f)

    task_args = TaskArguments(**pars)
    metadata_update = dummy_fail(**task_args.dict())

    if args.output:
        with open(args.output, "w") as fout:
            json.dump(metadata_update, fout)
    else:
        stdout.write(json.dumps(metadata_update))

    exit(0)
