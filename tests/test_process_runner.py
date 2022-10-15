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

from devtools import debug

import fractal_server.tasks.dummy as dummy
from fractal_server.app.models import Task
from fractal_server.app.models import Workflow
from fractal_server.app.runner import set_job_logger
from fractal_server.app.runner.process import process_workflow


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
        output_path=[tmp_path / "out.txt"],
        input_metadata={},
        logger=job_logger,
    )

    # check output
