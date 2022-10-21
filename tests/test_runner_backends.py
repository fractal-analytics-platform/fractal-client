"""
Copyright 2022 (C) Friedrich Miescher Institute for Biomedical Research and
University of Zurich

Original author(s):
Jacopo Nespolo <jacopo.nespolo@exact-lab.it>
Tommaso Comparin <tommaso.comparin@exact-lab.it>

This file is part of Fractal and was originally developed by eXact lab S.r.l.
<exact-lab.it> under contract with Liberali Lab from the Friedrich Miescher
Institute for Biomedical Research and Pelkmans Lab from the University of
Zurich.
"""
import logging

import pytest
from devtools import debug

import fractal_server.tasks.dummy as dummy
from fractal_server.app.models import Task
from fractal_server.app.models import Workflow
from fractal_server.app.runner import _backends
from fractal_server.app.runner import set_job_logger


@pytest.mark.parametrize(
    "backend",
    list(_backends.keys()),
)
async def test_runner(db, project_factory, MockCurrentUser, tmp_path, backend):
    """
    GIVEN a non-trivial workflow
    WHEN the workflow is processed
    THEN the tasks are correctly executed
    """
    debug(f"Testing with {backend=}")
    process_workflow = _backends[backend]

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
