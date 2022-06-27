from devtools import debug


PREFIX = "/api/v1/project"


async def test_project_creation(app, client, MockCurrentUser, db):
    debug(db)
    payload = dict(
        name="new project",
        project_dir="/some/path/",
    )
    res = await client.post(f"{PREFIX}/", json=payload)
    data = res.json()
    debug(data)
    assert res.status_code == 401

    with MockCurrentUser(sub="sub", scopes=["projects"]):
        debug(app.dependency_overrides)
        res = await client.post(f"{PREFIX}/", json=payload)
        data = res.json()
        debug(data)
        assert res.status_code == 201
