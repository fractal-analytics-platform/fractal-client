from devtools import debug


PREFIX = "/api/v1/project"


async def test_project_get(client, db, project_factory, MockCurrentUser):
    # unauthenticated

    res = await client.get(f"{PREFIX}/")
    assert res.status_code == 401

    # authenticated
    async with MockCurrentUser(persist=True) as user:
        other_project = await project_factory(user)

    async with MockCurrentUser(persist=True) as user:
        res = await client.get(f"{PREFIX}/")
        debug(res)
        assert res.status_code == 200
        assert res.json() == []

        await project_factory(user)
        res = await client.get(f"{PREFIX}/")
        data = res.json()
        debug(data)
        assert res.status_code == 200
        assert len(data) == 1

        project_id = data[0]["id"]
        res = await client.get(f"{PREFIX}/{project_id}")
        assert res.status_code == 200
        assert res.json()["id"] == project_id

        # fail on non existent project
        res = await client.get(f"{PREFIX}/666")
        assert res.status_code == 404

        # fail on other owner's project
        res = await client.get(f"{PREFIX}/{other_project.id}")
        assert res.status_code == 403


async def test_project_creation(app, client, MockCurrentUser, db):
    payload = dict(
        name="new project",
        project_dir="/some/path/",
    )
    res = await client.post(f"{PREFIX}/", json=payload)
    data = res.json()
    assert res.status_code == 401

    async with MockCurrentUser(persist=True):
        res = await client.post(f"{PREFIX}/", json=payload)
        data = res.json()
        assert res.status_code == 201
        debug(data)
        assert data["name"] == payload["name"]
        assert data["slug"] is not None
        assert data["project_dir"] == payload["project_dir"]
