import logging
from pathlib import Path
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
