import pytest


@pytest.mark.asyncio
async def test_alive(client):
    res = await client.get("/alive/")
    assert res.status_code == 200
    assert res.json()["alive"] is True
