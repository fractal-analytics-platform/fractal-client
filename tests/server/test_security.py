from devtools import debug


PREFIX = "/auth"


async def test_me(app, client, MockCurrentUser):
    # Anonymous
    res = await client.get(f"{PREFIX}/users/me")
    debug(res.json())
    assert res.status_code == 401
