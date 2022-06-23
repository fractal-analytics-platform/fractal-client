async def test_db_connection(db):
    # test active connection
    assert db.is_active
    # test bound
    assert db.get_bind()
