import pytest
from devtools import debug


@pytest.fixture
async def project_factory(db):
    from fractal.server.app.models import Project

    async def __project_factory(*kwargs):
        defaults = dict(name="project", project_dir="/tmp/")
        defaults.update(kwargs)
        project = Project(**defaults)
        db.add(project)
        await db.commit()
        await db.refresh(project)
        return project

    return __project_factory


PREFIX = "/api/v1/project"


async def test_project_get_list(client, db, project_factory):
    res = await client.get(f"{PREFIX}/")
    debug(res)
    assert res.json() == []

    await project_factory()
    res = await client.get(f"{PREFIX}/")
    data = res.json()
    debug(data)
    assert len(data) == 1


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
