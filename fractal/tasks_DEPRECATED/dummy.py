from pathlib import Path
from typing import Any
from typing import Dict
from typing import Iterable
from typing import Optional


def dummy(
    *,
    input_paths: Iterable[Path],
    output_path: Path,
    metadata: Optional[Dict[str, Any]] = None,
    component: Optional[Any] = None,
    # arguments of this task
    message: str,
    index: int = 0,
    **task_args,
) -> Path:
    """
    Dummy task

    This task appends to a json file the parameters it was called with, such
    that it is easy  to parse the file in a test settings.

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
        the arguments specific to the issue, i.e., `message` and `index` in the
        present exmaple

    **task_args (Any) :
        Task should always include a catch-all kwarg

    Retrun
    ------
    metadata_update (Dict[str, Any]) :
        a dictionary that will update the metadata
    """
    from datetime import datetime, timezone
    import json
    from json.decoder import JSONDecodeError

    if component:
        index = component

    payload = dict(
        task="DUMMY TASK",
        timestamp=datetime.now(timezone.utc).isoformat(),
        input_paths=[p.as_posix() for p in input_paths],
        output_path=output_path.as_posix(),
        metadata=metadata,
        message=message,
    )

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
    metadata_update = {"dummy": "dummy"}

    return metadata_update
