async def test_alive(client):
    res = await client.get("/alive/")
    data = res.json()
    assert res.status_code == 200
    assert data["alive"] is True
    assert data["deployment_type"] == "development"
