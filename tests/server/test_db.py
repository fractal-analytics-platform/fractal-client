from devtools import debug

from fractal.server.app.models.task import Task


async def test_db_connection(db):
    # test active connection
    assert db.is_active
    # test bound
    assert db.get_bind()

    from sqlmodel import select
    from fractal.server.app.models.security import UserOAuth

    stm = select(UserOAuth)
    res = await db.execute(stm)
    debug(res)


async def test_task_relations(db):
    t0 = Task(
        name="task0",
        import_name="task0",
        input_type="zarr",
        output_type="zarr",
    )
    t1 = Task(
        name="task1",
        import_name="task1",
        input_type="zarr",
        output_type="zarr",
    )
    db.add(t0)
    db.add(t1)
    t0.subtask_list.append(t1)
    await db.commit()
    await db.refresh(t0)
    await db.refresh(t1)
    debug(t0)
    debug(t1)
    assert len(t0.subtask_list) == 1
    assert t0.subtask_list[0] == t1
