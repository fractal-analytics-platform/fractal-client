import asyncio

from devtools import debug

from fractal.server.app.models import Subtask
from fractal.server.app.models import SubtaskRead
from fractal.server.app.models import TaskRead


async def test_task_relations(db, task_factory):
    """
    GIVEN two tasks
    WHEN a child task is added to a parent task
    THEN the relationship is correctly established in the database
    """
    parent = await task_factory(name="parent")
    child0 = await task_factory(name="child0")
    child1 = await task_factory(name="child1")

    db.add(parent)
    db.add(child0)
    child0_subtask = Subtask(parent=parent, subtask=child0)
    child1_subtask = Subtask(parent=parent, subtask=child1)
    db.add(child0_subtask)
    db.add(child1_subtask)
    await db.commit()
    await db.refresh(parent)
    await db.refresh(child0)
    await db.refresh(child0_subtask)
    await db.refresh(child1_subtask)
    debug(parent)
    assert child0_subtask.order == 0
    assert child1_subtask.order == 1
    assert len(parent.subtask_list) == 2
    assert parent.subtask_list[0].subtask == child0
    assert parent.subtask_list[0].subtask_id == child0.id

    subtask_read = SubtaskRead.from_orm(child0_subtask)
    child0_read = TaskRead.from_orm(child0)
    debug(subtask_read)
    assert subtask_read.subtask == child0_read

    # Swap children and check order
    child1_pop = parent.subtask_list.pop(1)
    parent.subtask_list.insert(0, child1_pop)
    db.add_all([parent, child0, child1])
    await db.commit()
    await asyncio.gather(
        *[db.refresh(item) for item in [parent, child0, child1]]
    )
    await db.commit()
    debug(parent)
    assert child0_subtask.order == 1
    assert child1_subtask.order == 0
