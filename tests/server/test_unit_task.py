from devtools import debug

from fractal.server.app.models import Subtask
from fractal.server.app.models import SubtaskRead
from fractal.server.app.models import Task
from fractal.server.app.models import TaskRead


async def test_task_relations(db, task_factory):
    """
    GIVEN two tasks
    WHEN a child task is added to a parent task
    THEN the relationship is correctly established in the database
    """
    parent = Task(
        name="parent",
        resource_type="task",
        module="parent",
        input_type="zarr",
        output_type="zarr",
    )
    child = Task(
        name="child",
        resource_type="task",
        module="child",
        input_type="zarr",
        output_type="zarr",
    )
    db.add(parent)
    db.add(child)
    child_subtask = Subtask(parent=parent, subtask=child)
    db.add(child_subtask)
    await db.commit()
    await db.refresh(parent)
    await db.refresh(child)
    await db.refresh(child_subtask)
    debug(parent)
    debug(child)
    assert len(parent.subtask_list) == 1
    assert parent.subtask_list[0].subtask == child
    assert parent.subtask_list[0].subtask_id == child.id

    subtask_read = SubtaskRead.from_orm(child_subtask)
    child_read = TaskRead.from_orm(child)
    debug(subtask_read)
    assert subtask_read.subtask == child_read
