import pytest
from devtools import debug


@pytest.fixture
async def project_factory(db):
    from fractal.server.app.models import Project

    async def __project_factory(user, **kwargs):
        defaults = dict(
            name="project",
            project_dir="/tmp/",
            user_owner_id=user.id,
            slug="slug",
        )
        defaults.update(kwargs)
        project = Project(**defaults)
        db.add(project)
        await db.commit()
        await db.refresh(project)
        return project

    return __project_factory


PREFIX = "/api/v1/project"


async def test_project_get_list(client, db, project_factory, MockCurrentUser):
    # unauthenticated

    res = await client.get(f"{PREFIX}/")
    assert res.status_code == 401

    # authenticated
    async with MockCurrentUser(persist=True) as user:
        res = await client.get(f"{PREFIX}/")
        debug(res)
        assert res.json() == []

        await project_factory(user)
        res = await client.get(f"{PREFIX}/")
        data = res.json()
        debug(data)
        assert len(data) == 1


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
