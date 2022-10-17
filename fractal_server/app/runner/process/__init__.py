import logging
import subprocess  # nosec
from pathlib import Path
from shlex import split as shlex_split
from typing import Any
from typing import Dict
from typing import List

from ...models import Workflow

"""
Process Bakend

This backend runs fractal workflows as separate processes using a python
process pool.

Incidentally, it represents the reference implementation for a backend.
"""


def _call_command_wrapper(cmd: str):
    """
    Call command and return stdout, stderr, retcode
    """

    return subprocess.run(shlex_split(cmd), capture_output=True)  # nosec


async def process_workflow(
    *,
    workflow: Workflow,
    input_paths: List[Path],
    output_path: Path,
    input_metadata: Dict[str, Any],
    logger: logging.Logger,
    username: str = None,
) -> Dict[str, Any]:
    raise NotImplementedError
