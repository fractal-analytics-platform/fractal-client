from devtools import debug


async def test_db_connection(db):
    # test active connection
    assert db.is_active
    # test bound
    assert db.get_bind()
    debug(db.get_bind())
    debug(db.get_bind().url.database)
    assert db.get_bind().url.database is not None

    from sqlmodel import select
    from fractal_server.app.models.security import UserOAuth

    stm = select(UserOAuth)
    res = await db.execute(stm)
    debug(res)


async def test_sync_db(db_sync, db):
    """
    GIVEN a database and a sync and an async connections to it
    WHEN crud operations are executed with one connection
    THEN results are consistent with the other connection
    """
    assert db_sync.is_active
    assert db_sync.get_bind()

    from sqlmodel import select
    from fractal_server.app.models.task import Task

    db.add(
        Task(
            name="mytask",
            input_type="image",
            output_type="zarr",
            command="cmd",
            source="/source",
        )
    )
    await db.commit()

    # Async
    stm = select(Task)
    res = await db.execute(stm)
    task_list = res.scalars().all()
    assert len(task_list) == 1
    assert task_list[0].name == "mytask"

    # Sync
    res = db_sync.execute(stm)
    task_list = res.scalars().all()
    assert len(task_list) == 1
    assert task_list[0].name == "mytask"
