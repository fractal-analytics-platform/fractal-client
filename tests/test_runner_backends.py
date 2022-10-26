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

from fractal_server.app.models import Task
from fractal_server.app.models import Workflow
from fractal_server.app.runner import _backends
from fractal_server.app.runner.common import close_job_logger
from fractal_server.app.runner.common import set_job_logger
from fractal_server.tasks import dummy
from fractal_server.tasks import dummy_parallel


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

    tp = Task(
        name="dummy_parallel",
        command=f"python {dummy_parallel.__file__}",
        source=dummy.__file__,
        input_type="Any",
        output_type="Any",
        default_args=dict(parallelization_level="index"),
    )

    # Create a workflow with the dummy task as member
    wf = Workflow(name="wf", project_id=prj.id)

    db.add_all([tk, tp, wf])
    await db.commit()
    await db.refresh(tk)
    await db.refresh(tp)
    await db.refresh(wf)

    await wf.insert_task(tk.id, db=db, args=dict(message="task 0"))
    await wf.insert_task(tk.id, db=db, args=dict(message="task 1"))
    await wf.insert_task(tp.id, db=db, args=dict(message="task 2"))
    await db.refresh(wf)

    debug(tk)
    debug(wf)

    # process workflow
    logger_name = "job_logger"
    logger = set_job_logger(
        logger_name=logger_name,
        log_file_path=tmp_path / "job.log",
        level=logging.DEBUG,
    )
    out = await process_workflow(
        workflow=wf,
        input_paths=[tmp_path / "*.txt"],
        output_path=tmp_path / "out.json",
        input_metadata={},
        logger_name=logger_name,
        workflow_dir=tmp_path,
    )
    close_job_logger(logger)
    debug(out)
    assert "dummy" in out.metadata
    assert "dummy" in out.metadata
    assert out.metadata["history"] == [
        tk.name,
        tk.name,
        f"{tp.name}: [0, 1, 2]",
    ]
