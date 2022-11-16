"""
Copyright 2022 (C) Friedrich Miescher Institute for Biomedical Research and
University of Zurich

Original author(s):
Jacopo Nespolo <jacopo.nespolo@exact-lab.it>

This file is part of Fractal and was originally developed by eXact lab S.r.l.
<exact-lab.it> under contract with Liberali Lab from the Friedrich Miescher
Institute for Biomedical Research and Pelkmans Lab from the University of
Zurich.
"""
from os import environ

import pytest
from devtools import debug

PREFIX = "/api/v1"

environ["RUNNER_MONITORING"] = "0"


@pytest.fixture
async def collect_tasks(MockCurrentUser, client, dummy_task_package):
    async with MockCurrentUser(persist=True):
        # COLLECT DUMMY TASKS
        res = await client.post(
            f"{PREFIX}/task/collect/pip/",
            json=dict(package=dummy_task_package.as_posix()),
        )
        state = res.json()
        data = state["data"]

        res = await client.get(f"{PREFIX}/task/collect/{state['id']}")
        assert res.status_code == 200
        state = res.json()
        data = state["data"]
        task_list = data["task_list"]
    return task_list


async def test_full_workflow(
    client,
    MockCurrentUser,
    testdata_path,
    tmp_path,
    collect_tasks,
    project_factory,
    dataset_factory,
):

    async with MockCurrentUser(persist=True) as user:
        project = await project_factory(user)
        debug(project)
        project_id = project.id
        input_dataset = await dataset_factory(
            project, name="input", type="image", read_only=True
        )
        input_dataset_id = input_dataset.id

        # EDIT DEFAULT DATASET TO SET TYPE IMAGE

        res = await client.patch(
            f"{PREFIX}/project/{project_id}/{input_dataset_id}",
            json={"type": "image", "read_only": True},
        )
        debug(res.json())
        assert res.status_code == 200

        # ADD TEST IMAGES AS RESOURCE TO INPUT DATASET

        res = await client.post(
            f"{PREFIX}/project/{project_id}/{input_dataset_id}",
            json={
                "path": (testdata_path / "png").as_posix(),
                "glob_pattern": "*.png",
            },
        )
        debug(res.json())
        assert res.status_code == 201

        # CREATE OUTPUT DATASET AND RESOURCE

        res = await client.post(
            f"{PREFIX}/project/{project_id}/",
            json=dict(
                name="output dataset",
                type="json",
            ),
        )
        debug(res.json())
        assert res.status_code == 201
        output_dataset = res.json()
        output_dataset_id = output_dataset["id"]

        res = await client.post(
            f"{PREFIX}/project/{project_id}/{output_dataset['id']}",
            json=dict(path=tmp_path.as_posix(), glob_pattern="out.json"),
        )
        out_resource = res.json()
        debug(out_resource)
        assert res.status_code == 201

        # CHECK WHERE WE ARE AT
        res = await client.get(f"{PREFIX}/project/{project_id}")
        debug(res.json())
        project_dict = res.json()

        # CREATE WORKFLOW
        res = await client.post(
            f"{PREFIX}/workflow/",
            json=dict(name="test workflow", project_id=project_dict["id"]),
        )
        debug(res.json())
        assert res.status_code == 201
        workflow_dict = res.json()
        workflow_id = workflow_dict["id"]

        # Add a dummy task
        res = await client.post(
            f"{PREFIX}/workflow/{workflow_id}/add-task/",
            json=dict(task_id=collect_tasks[0]["id"]),
        )
        debug(res.json())
        assert res.status_code == 201

        # Add another (identical) dummy task, to make sure that this is allowed
        res = await client.post(
            f"{PREFIX}/workflow/{workflow_id}/add-task/",
            json=dict(task_id=collect_tasks[0]["id"]),
        )
        debug(res.json())
        assert res.status_code == 201

        # Add a dummy_parallel task
        res = await client.post(
            f"{PREFIX}/workflow/{workflow_id}/add-task/",
            json=dict(task_id=collect_tasks[1]["id"]),
        )
        debug(res.json())
        assert res.status_code == 201

        # EXECUTE WORKFLOW

        payload = dict(
            project_id=project_id,
            input_dataset_id=input_dataset_id,
            output_dataset_id=output_dataset_id,
            workflow_id=workflow_id,
            overwrite_input=False,
        )
        debug(payload)
        res = await client.post(
            f"{PREFIX}/project/apply/",
            json=payload,
        )
        debug(res.json())
        assert res.status_code == 202

        # Verify output
        res = await client.get(
            f"{PREFIX}/dataset/{project_id}/{output_dataset_id}"
        )
        data = res.json()
        debug(data)
        assert "dummy" in data["meta"]
