"""
Copyright 2022 (C) Friedrich Miescher Institute for Biomedical Research and
University of Zurich

Original author(s):
Jacopo Nespolo <jacopo.nespolo@exact-lab.it>

This file is part of Fractal and was originally developed by eXact lab S.r.l.
<exact-lab.it> under contract with Liberali Lab from the Friedrich Miescher
Institute for Biomedical Research and Pelkmans Lab from the University of
Zurich.
"""
import json
import logging
from concurrent.futures import ThreadPoolExecutor
from typing import Dict

from devtools import debug
from pydantic import BaseModel

import fractal_server.tasks.dummy as dummy
from fractal_server.app.models import Task
from fractal_server.app.models import Workflow
from fractal_server.app.runner import set_job_logger
from fractal_server.app.runner.common import TaskParameters
from fractal_server.app.runner.process import _call_command_wrapper
from fractal_server.app.runner.process import call_single_task
from fractal_server.app.runner.process import process_workflow
from fractal_server.app.runner.process import recursive_task_submission
from fractal_server.tasks import dummy as dummy_module
from fractal_server.tasks import dummy_parallel as dummy_parallel_module


class MockTask(BaseModel):
    command: str
    parallelization_level: str = None


class MockWorkflowTask(BaseModel):
    order: int = 0
    task: MockTask
    arguments: Dict = {}


MOCKPARALLELTASK_NAME = "This is just a name"


class MockParallelTask(BaseModel):
    name: str = MOCKPARALLELTASK_NAME
    command: str
    parallelization_level: str


class MockParallelWorkflowTask(BaseModel):
    order: int = 0
    task: MockParallelTask
    arguments: Dict = {}


async def test_command_wrapper():
    with ThreadPoolExecutor() as executor:
        future = executor.submit(
            _call_command_wrapper, f"ls -al {dummy.__file__}"
        )
    result = future.result()
    debug(result.stdout, result.stderr, result.returncode)
    assert dummy.__file__ in result.stdout.decode("utf-8")


def test_call_single_task(tmp_path):
    task = MockWorkflowTask(
        task=MockTask(command=f"python {dummy_module.__file__}"),
        arguments=dict(message="test"),
        order=0,
    )
    task_pars = TaskParameters(
        input_paths=[tmp_path],
        output_path=tmp_path,
        metadata={},
        logger=logging.getLogger(),
    )

    debug(task)

    out = call_single_task(
        task=task, task_pars=task_pars, workflow_dir=tmp_path
    )
    debug(out)
    assert isinstance(out, TaskParameters)
    # check specific of the dummy task
    assert out.metadata["dummy"] == "dummy 0"


def test_recursive_task_submission_step0(tmp_path):
    """
    GIVEN a workflow with a single task
    WHEN it is passed to the recursive task submission
    THEN it is correctly executed, i.e., step 0 of the induction
    """
    INDEX = 666
    task_list = [
        MockWorkflowTask(
            task=MockTask(command=f"python {dummy_module.__file__}"),
            arguments=dict(message="test", index=INDEX),
            order=0,
        )
    ]
    job_logger = set_job_logger(
        logger_name="job_logger",
        log_file_path=tmp_path / "job.log",
        level=logging.DEBUG,
    )
    task_pars = TaskParameters(
        input_paths=[tmp_path],
        output_path=tmp_path,
        metadata={},
        logger=job_logger,
    )

    with ThreadPoolExecutor() as executor:
        res = recursive_task_submission(
            executor=executor,
            task_list=task_list,
            task_pars=task_pars,
            workflow_dir=tmp_path,
        )
        debug(res.result())
        assert res.result().metadata["dummy"] == f"dummy {INDEX}"


def test_recursive_parallel_task_submission_step0(tmp_path):
    """
    GIVEN a workflow with a single parallel task
    WHEN it is passed to the recursive task submission
    THEN it is correctly executed, i.e., step 0 of the induction
    """
    LIST_INDICES = ["0", "1"]
    MESSAGE = "test message"
    task_list = [
        MockParallelWorkflowTask(
            task=MockParallelTask(
                command=f"python {dummy_parallel_module.__file__}",
                parallelization_level="index",
            ),
            arguments=dict(message=MESSAGE),
            order=0,
        )
    ]
    job_logger = set_job_logger(
        logger_name="job_logger",
        log_file_path=tmp_path / "job.log",
        level=logging.DEBUG,
    )
    output_path = tmp_path / "output/*.json"
    task_pars = TaskParameters(
        input_paths=[tmp_path],
        output_path=output_path,
        metadata={"index": LIST_INDICES},
        logger=job_logger,
    )

    debug(task_list)
    debug(task_pars)

    with ThreadPoolExecutor() as executor:
        res = recursive_task_submission(
            executor=executor,
            task_list=task_list,
            task_pars=task_pars,
            workflow_dir=tmp_path,
        )
        debug(res.result())
        assert MOCKPARALLELTASK_NAME in res.result().metadata["history"][0]

    # Validate results
    assert output_path.parent.exists()
    output_files = list(output_path.parent.glob("*"))
    debug(output_files)
    assert len(output_files) == len(LIST_INDICES)

    for output_file in output_files:
        with output_file.open("r") as fin:
            data = json.load(fin)
        assert output_file.name == f'{data["component"]}.json'
        assert data["message"] == MESSAGE


def test_recursive_parallel_task_submission_inductive_step(tmp_path):
    """
    GIVEN a workflow with three global/parallel/global tasks
    WHEN it is passed to the recursive task submission
    THEN it is correctly executed, i.e., n => n+1
    """
    # FIXME
    pass


def test_recursive_task_submission_inductive_step(tmp_path):
    """
    GIVEN a workflow with two or more tasks
    WHEN it is passed to the recursive task submission
    THEN it is correctly executed, i.e., n => n+1
    """
    METADATA_0 = {}
    METADATA_1 = dict(dummy="dummy 0")  # for how dummy task works

    task_list = [
        MockWorkflowTask(
            task=MockTask(command=f"python {dummy_module.__file__}"),
            arguments=dict(message="test 0", index=0),
            order=0,
        ),
        MockWorkflowTask(
            task=MockTask(command=f"python {dummy_module.__file__}"),
            arguments=dict(message="test 1", index=1),
            order=1,
        ),
    ]
    job_logger = set_job_logger(
        logger_name="job_logger",
        log_file_path=tmp_path / "job.log",
        level=logging.DEBUG,
    )
    task_pars = TaskParameters(
        input_paths=[tmp_path],
        output_path=tmp_path / "output.json",
        metadata=METADATA_0,
        logger=job_logger,
    )

    with ThreadPoolExecutor() as executor:
        res = recursive_task_submission(
            executor=executor,
            task_list=task_list,
            task_pars=task_pars,
            workflow_dir=tmp_path,
        )

    output = res.result()
    debug(output)
    with open(output.output_path, "r") as f:
        data = json.load(f)
    debug(data)
    assert len(data) == 2
    assert data[0]["metadata"] == METADATA_0
    assert data[1]["metadata"] == METADATA_1


async def test_runner(db, project_factory, MockCurrentUser, tmp_path):
    """
    GIVEN a non-trivial workflow
    WHEN the workflow is processed
    THEN the tasks are correctly executed
    """
    async with MockCurrentUser(persist=True) as user:
        prj = await project_factory(user=user)

    # Add dummy task as a Task
    tk = Task(
        name="dummy",
        command=f"python {dummy.__file__}",
        source=dummy.__file__,
        input_type="Any",
        output_type="Any",
    )

    # Create a workflow with the dummy task as member
    wf = Workflow(name="wf", project_id=prj.id)

    db.add_all([tk, wf])
    await db.commit()
    await db.refresh(tk)
    await db.refresh(wf)

    await wf.insert_task(tk, db=db, args=dict(message="task 0"))
    await wf.insert_task(tk, db=db, args=dict(message="task 1"))
    await db.refresh(wf)

    debug(tk)
    debug(wf)

    # process workflow
    job_logger = set_job_logger(
        logger_name="job_logger",
        log_file_path=tmp_path / "job.log",
        level=logging.DEBUG,
    )
    out = await process_workflow(
        workflow=wf,
        input_paths=[tmp_path / "*.txt"],
        output_path=tmp_path / "out.json",
        input_metadata={},
        logger=job_logger,
        workflow_dir=tmp_path,
    )
    debug(out)
    assert "dummy" in out.metadata
