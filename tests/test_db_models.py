from sqlmodel import select

from fractal_server.app.models import Project


async def test_project_name_not_unique(MockCurrentUser, db, project_factory):
    """
    GIVEN the fractal_server database
    WHEN I create two projects with the same name and same user
    THEN no exception is raised
    """
    PROJ_NAME = "project name"
    async with MockCurrentUser(persist=True) as user:
        p0 = await project_factory(user, name=PROJ_NAME)
        p1 = await project_factory(user, name=PROJ_NAME)

    stm = select(Project).where(Project.name == PROJ_NAME)
    res = await db.execute(stm)
    project_list = res.scalars().all()
    assert len(project_list) == 2
    assert p0 in project_list
    assert p1 in project_list
