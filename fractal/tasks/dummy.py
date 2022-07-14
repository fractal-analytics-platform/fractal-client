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
    # arguments of this task
    message: str,
    **task_args,
):
    """
    Dummy task

    This task appends to a json file the parameters it was called with, such
    that it is easy  to parse the file in a test settings.

    Incidentally, this task defines the reference implementation of a task.
    """
    from datetime import datetime, timezone
    import json
    from json.decoder import JSONDecodeError

    payload = dict(
        task="DUMMY TASK",
        timestamp=datetime.now(timezone.utc).isoformat(),
        input_path=[p.as_posix() for p in input_paths],
        output_path=output_path.as_posix(),
        metadata=metadata,
        message=message,
    )

    with open(output_path, "r+") as fout:
        try:
            data = json.load(fout)
        except JSONDecodeError:
            data = []
        data.append(payload)
        json.dump(data, fout)
