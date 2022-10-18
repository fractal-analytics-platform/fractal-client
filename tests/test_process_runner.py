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
import logging
from concurrent.futures import ThreadPoolExecutor

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
from fractal_server.tasks import dummy as dummy_module


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

    await wf.insert_task(tk, db=db)
    await db.refresh(wf)

    debug(tk)
    debug(wf)

    # process workflow
    job_logger = set_job_logger(
        logger_name="job_logger",
        log_file_path=tmp_path / "job.log",
        level=logging.INFO,
    )
    await process_workflow(
        workflow=wf,
        input_paths=[tmp_path / "*.txt"],
        output_path=tmp_path / "out.txt",
        input_metadata={},
        logger=job_logger,
    )

    # check output
    raise NotImplementedError


async def test_command_wrapper():
    with ThreadPoolExecutor() as executor:
        future = executor.submit(
            _call_command_wrapper, f"ls -al {dummy.__file__}"
        )
    result = future.result()
    debug(result.stdout, result.stderr, result.returncode)
    assert dummy.__file__ in result.stdout.decode("utf-8")


def test_recursive_task_submission():
    from pydantic import BaseModel

    class TaskParameters(BaseModel):
        input_path: str
        output_path: str
        metadata: str

    def recursive_task_submission(task_list, task_parameters):
        try:
            *depenency_tasks, this_task = task_list
        except ValueError:
            return task_parameters
        dependency = recursive_task_submission(
            depenency_tasks, task_parameters
        )

        this_task_parameters = dependency
        debug(this_task, this_task_parameters)
        return TaskParameters(
            input_path=task_parameters.output_path,
            output_path=task_parameters.output_path,
            metadata=f"metaupdate {this_task}",
        )

    recursive_task_submission(
        ["task0", "task1", "task2"],
        task_parameters=TaskParameters(
            input_path="input0",
            output_path="output",
            metadata="metadata0",
        ),
    )
    assert False


def test_call_single_task(tmp_path):
    class MockTask(BaseModel):
        command: str

    class MockWorkflowTask(BaseModel):
        order: int = 0
        task: MockTask
        arguments = {}

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
    debug(out.stderr.decode("utf-8"))
    debug(out.returncode)
    assert False
