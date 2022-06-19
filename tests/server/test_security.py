from devtools import debug


async def test_token_endpoint(client):
    res = await client.post(
        "/v1/token", data=dict(username="username", password="password")
    )
    debug(res.json())
    assert res.status_code == 200

    res = await client.post(
        "/v1/token", data=dict(username="fail", password="fail")
    )
    debug(res.json())
    assert res.status_code == 401


async def test_who_am_i(client):
    # Anonymous
    res = await client.get("v1/me")
    assert res.status_code == 401

    res = await client.post(
        "/v1/token", data=dict(username="username", password="password")
    )
    data = res.json()
    headers = dict(
        Authorization=f"{data['token_type']} {data['access_token']}"
    )
    res = await client.get("v1/me", headers=headers)
    data = res.json()
    debug(data)
    assert res.status_code == 200
    assert data["sub"] == "1234"
    assert data["name"] == "User Name"
