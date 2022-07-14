import json
from pathlib import Path
from tempfile import NamedTemporaryFile

import pytest
from devtools import debug
from sqlmodel import select

from fractal.server.app.models import Subtask
from fractal.server.app.models import Task
from fractal.server.app.models import TaskRead
from fractal.server.app.runner import _atomic_task_factory
from fractal.server.app.runner import submit_workflow


DUMMY_MESSAGE = "dummy task message"
DUMMY_SUBTASK_MESSAGE = "dummy subtask message"

dummy_task = Task(
    name="dummy",
    resource_type="core task",
    module="fractal.tasks.dummy:dummy",
    default_args={"message": DUMMY_MESSAGE},
)
dummy_subtask = Subtask(
    subtask=dummy_task, args={"message": DUMMY_SUBTASK_MESSAGE}
)


@pytest.mark.parametrize(
    ("task", "message"),
    [(dummy_task, DUMMY_MESSAGE), (dummy_subtask, DUMMY_SUBTASK_MESSAGE)],
)
async def test_atomic_task_factory(task, message):
    """
    GIVEN
        * a task or subtask
        * input_path iterable
        * output_path
    WHEN passed to the task factory
    THEN
        * the relative PARSL workflow is correctly
        * it can run
        * the output is as expected
    """
    input_path_str = "/input/path"
    output_file = NamedTemporaryFile()

    parsl_app = _atomic_task_factory(
        task=task,
        input_paths=[Path(input_path_str)],
        output_path=Path(output_file.name),
        metadata={},
    )

    debug(parsl_app)
    parsl_app.result()

    data = json.load(output_file)
    debug(data)
    assert len(data) == 1
    assert data[0]["message"] == message


async def test_apply_workflow(
    db,
    client,
    collect_tasks,
    MockCurrentUser,
    project_factory,
    dataset_factory,
    resource_factory,
    task_factory,
):
    """
    GIVEN
        * an input dataset and relative resource(s)
        * an output dataset and relative resource
        * a non-trivial workflow
    WHEN one applys the workflow to the input dataset
    THEN
        * the workflow is executed correctly
        * the output is correctly written in the output resource
    """

    # CREATE RESOURCES
    async with MockCurrentUser(persist=True) as user:
        prj = await project_factory(user)
        ds = await dataset_factory(prj, type="image")
        out_ds = await dataset_factory(prj, type="image", name="out_ds")

        resource = await resource_factory(ds)

        debug(ds)
        debug(resource)

    # CREATE NONTRIVIAL WORKFLOW
    wf = await task_factory(
        name="worfklow",
        module=None,
        resource_type="workflow",
        input_type="image",
    )
    debug(wf)

    stm = select(Task).where(Task.name == "dummy")
    res = await db.execute(stm)
    dummy_task = res.scalar()

    await wf.add_subtask(db, subtask=dummy_task)
    debug(TaskRead.from_orm(wf))

    await submit_workflow(input_dataset=ds, output_dataset=out_ds, workflow=wf)
