import pytest
from devtools import debug
from sqlmodel import select

from fractal_server.app.models import Task
from fractal_server.app.models import TaskCreate
from fractal_server.tasks import collect_tasks

async def test_workflow_post(db, client, MockCurrentUser):
    fake_workflow = {
        "name": "My Workflow",
        "resource_type": "core task",
        "input_type": "Boh",
        "output_type": "Boh",
    }
    async with MockCurrentUser(persist=True):
        
        res1 = await client.post(
            "api/v1/workflow/",
            json=fake_workflow,
        )
        res2 = await client.post(
            "api/v1/workflow/",
            json=fake_workflow,
        )
        assert res1.status_code==201
        assert res2.status_code==422
    
    debug(res1, res1.json())
    debug(res2, res2.json())
        