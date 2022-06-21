from devtools import debug


PREFIX = "/v1/project"


async def test_project_creation(app, client, MockCurrentUser):
    payload = dict(
        name="new project",
        project_dir="/some/path/",
    )
    res = await client.post(f"{PREFIX}/", json=payload)
    data = res.json()
    debug(data)
    assert res.status_code == 201
