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
import pytest
from devtools import debug


PREFIX = "/api/v1"


@pytest.mark.xfail
async def test_project_creation(
    app, client, MockCurrentUser, db, testdata_path, collect_tasks
):
    async with MockCurrentUser(persist=True):

        # CREATE PROJECT

        res = await client.post(
            f"{PREFIX}/project/",
            json=dict(
                name="test project",
                project_dir="/tmp/",
            ),
        )
        assert res.status_code == 201
        project = res.json()
        project_id = project["id"]
        dataset_id = project["dataset_list"][0]["id"]

        # ADD RESOURCE TO DATASET

        res = await client.post(
            f"{PREFIX}/project/{project_id}/{dataset_id}",
            json={"path": testdata_path.as_posix()},
        )
        debug(res.json())
        assert res.status_code == 201

        # GET WORKFLOW ID (in this case just a dummy task)
        res = await client.get(f"{PREFIX}/task/")
        assert res.status_code == 200
        data = res.json()
        debug(data)
        task = data[0]
        workflow_id = task["id"]

        # EXECUTE WORKFLOW

        payload = dict(
            project_id=project_id,
            workflow_id=workflow_id,
            input_dataset_id=dataset_id,
        )
        res = await client.post(
            f"{PREFIX}/project/apply/",
            json=payload,
        )
        debug(res.json())
        assert res.status_code == 202
