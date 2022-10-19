from re import A
import pytest
from devtools import debug
from sqlmodel import select

from fractal_server.app.models import Workflow
from fractal_server.app.models import WorkflowRead
from fractal_server.app.models import WorkflowCreate
from fractal_server.app.models import LinkTaskWorkflow

from fractal_server.app.models import Project
from fractal_server.app.models import Task

async def test_workflow_post(
    db, client, MockCurrentUser, project_factory, task_factory
):
    async with MockCurrentUser(persist=True) as user:
        res = await client.post(
            "api/v1/workflow/",
            json={"name":"My Workflow"}, 
        )
        assert res.status_code==422 # no project_id

        project = await project_factory(user)
        p_id = project.id
        workflow = WorkflowCreate(
            name="My Workflow",
            project_id=p_id,
        )

        res = await client.post(
            "api/v1/workflow/",
            json=workflow.dict(),
        )
        assert res.status_code==201

        stm = (
            select(Workflow)
            .where(Workflow.project_id == p_id)
        )
        _workflow = await db.execute(stm)
        workflow = _workflow.scalars().one()

        t0 = await task_factory()
        workflow.insert_task(t0)
        t1 = await task_factory()
        workflow.insert_task(t1)

        db.merge(workflow)
        db.commit()
        db.refresh(workflow)
        w_id = workflow.id
        
        stm = select(Project)
        _project = await db.execute(stm)
        project_list = _project.scalars().all()
        assert len(project_list)==1
        project = project_list[0]

        stm = select(Task)
        _task = await db.execute(stm)
        task_list = _task.scalars().all()
        assert len(task_list)==2
        task_id_list = [task.id for task in task_list]

        stm = select(LinkTaskWorkflow)
        _link = await db.execute(stm)
        link_list = _link.scalars().all()
        
        for link in link_list:
            assert link.workflow_id==w_id
            assert link.task_id in task_id_list
            task_id_list.remove(link.task_id)
        assert task_id_list == []
