from uuid import uuid4

import pytest
from fastapi import HTTPException

from fractal_server.app.api.v1.project import get_project_check_owner


async def test_proejct_membership(db, project_factory, MockCurrentUser):
    async with MockCurrentUser(persist=True) as user:
        prj = await project_factory(user=user)

        # existing project and member user
        out_prj = await get_project_check_owner(
            project_id=prj.id, user_id=user.id, db=db
        )
        assert out_prj == prj

        # non-existing project
        with pytest.raises(HTTPException) as e:
            out_prj = await get_project_check_owner(
                project_id=666, user_id=user.id, db=db
            )
        assert e.value.status_code == 404

        # non-member user
        with pytest.raises(HTTPException) as e:
            out_prj = await get_project_check_owner(
                project_id=prj.id, user_id=uuid4(), db=db
            )
        assert e.value.status_code == 403
