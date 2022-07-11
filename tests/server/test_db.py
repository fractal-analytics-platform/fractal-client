from devtools import debug


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
