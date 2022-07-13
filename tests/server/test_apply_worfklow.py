from devtools import debug
from sqlmodel import select

from fractal.server.app.models import Task


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
    GIVEN a dataset and a non-trivial workflow
    WHEN one applys the workflow to the dataset
    THEN the job is correctly submitted and executed
    """

    # CREATE RESOURCES
    async with MockCurrentUser(persist=True) as user:
        prj = await project_factory(user)
        ds = await dataset_factory(prj, type="image")

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

    debug(dummy_task)
