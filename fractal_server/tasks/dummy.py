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
from datetime import datetime
from datetime import timezone
from json.decoder import JSONDecodeError
from pathlib import Path
from sys import stdout
from typing import Any
from typing import Dict
from typing import Iterable
from typing import List
from typing import Optional

from pydantic import BaseModel


logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)


def dummy(
    *,
    input_paths: Iterable[Path],
    output_path: Path,
    metadata: Optional[Dict[str, Any]] = None,
    # arguments of this task
    message: str,
    index: int = 0,
    raise_error: bool = False,
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

    Retrun
    ------
    metadata_update (Dict[str, Any]) :
        a dictionary that will update the metadata
    """
    logger.info("ENTERING dummy task")

    if raise_error:
        raise ValueError(message)

    payload = dict(
        task="DUMMY TASK",
        timestamp=datetime.now(timezone.utc).isoformat(),
        input_paths=[p.as_posix() for p in input_paths],
        output_path=output_path.as_posix(),
        metadata=metadata,
        message=message,
    )

    if not output_path.parent.is_dir():
        output_path.parent.mkdir()

    if not output_path.as_posix().endswith(".json"):
        filename_out = f"{index}.json"
        out_fullpath = output_path / filename_out
    else:
        out_fullpath = output_path

    try:
        with open(out_fullpath, "r") as fin:
            data = json.load(fin)
    except (JSONDecodeError, FileNotFoundError):
        data = []
    data.append(payload)
    with open(out_fullpath, "w") as fout:
        json.dump(data, fout, indent=2)

    # Update metadata
    metadata_update = {"dummy": f"dummy {index}"}

    logger.info("EXITING dummy task")

    return metadata_update


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
        message: str
        index: int = 0
        raise_error: bool = False

    parser = ArgumentParser()
    parser.add_argument("-j", "--json", help="Read parameters from json file")
    parser.add_argument(
        "--metadata-out",
        help=(
            "Output file to redirect serialised returned data "
            "(default stdout)"
        ),
    )

    args = parser.parse_args()

    if args.metadata_out and Path(args.metadata_out).exists():
        logger.error(
            f"Output file {args.metadata_out} already exists. Terminating"
        )
        exit(1)

    pars = {}
    if args.json:
        with open(args.json, "r") as f:
            pars = json.load(f)

    task_args = TaskArguments(**pars)
    metadata_update = dummy(**task_args.dict())

    if args.metadata_out:
        with open(args.metadata_out, "w") as fout:
            json.dump(metadata_update, fout)
    else:
        stdout.write(json.dumps(metadata_update))

    exit(0)
