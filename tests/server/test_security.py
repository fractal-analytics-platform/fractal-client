from devtools import debug


async def test_token_endpoint(client):
    res = await client.post(
        "/v1/token", data=dict(username="user", password="password")
    )
    debug(res.json())
    assert res.status_code == 200
