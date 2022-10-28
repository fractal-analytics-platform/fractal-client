import json
import logging
from concurrent.futures import Future

import pytest
from devtools import debug

from .fixtures_tasks import MockTask
from .fixtures_tasks import MockWorkflowTask
from fractal_server.app.runner.common import close_job_logger
from fractal_server.app.runner.common import TaskParameters
from fractal_server.tasks import dummy as dummy_module
from fractal_server.utils import set_logger


try:
    import parsl  # noqa: F401
    from fractal_server.app.runner._parsl import _serial_task_assembly
    from fractal_server.app.runner._parsl._setup import load_parsl_config
except ImportError:
    pytest.skip(allow_module_level=True)


def test_import_parsl_backend(unset_deployment_type):
    """
    Test that we can import _parsl even though the settings would fail the
    check
    """
    import fractal_server.app.runner._parsl  # noqa: F401


def test_unit_serial_task_assembly(tmp_path):

    INDEX = 666
    MESSAGE = "this message"
    EXECUTORS = "all"

    workflow_task = MockWorkflowTask(
        task=MockTask(name="task0", command=f"python {dummy_module.__file__}"),
        arguments=dict(message=MESSAGE, index=INDEX),
        order=0,
    )

    logger_name = "test_logger"
    logger = set_logger(
        logger_name=logger_name,
        log_file_path=tmp_path / "task.log",
        level=logging.DEBUG,
    )
    task_pars = TaskParameters(
        input_paths=[tmp_path],
        output_path=tmp_path,
        metadata={},
        logger=logging.getLogger(),
    )
    task_pars_depend_future = Future()
    task_pars_depend_future.set_result(task_pars)

    with load_parsl_config(
        workflow_id=42,  # here
        workflow_name="workflow_name",  # here
        workflow_dir=tmp_path,
        username=None,
        logger_name=logger_name,
    ) as dfk:
        out = _serial_task_assembly(
            data_flow_kernel=dfk,
            task=workflow_task,
            task_pars_depend_future=task_pars_depend_future,
            workflow_dir=tmp_path,
            executors=EXECUTORS,
        )
        debug(out.result())
    close_job_logger(logger)

    # metadata file exists
    output_file = tmp_path / "0.metadiff.json"
    with output_file.open("r") as fout:
        data = json.load(fout)

    assert data["dummy"] == f"dummy {INDEX}"
    # output file exists
    output_file = tmp_path / f"{INDEX}.json"
    with output_file.open("r") as fout:
        data = json.load(fout)

    assert len(data) == 1
    assert data[0]["message"] == MESSAGE
