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


PREFIX = "/api/v1"


async def test_project_creation(
    app, client, MockCurrentUser, db, testdata_path
):
    with MockCurrentUser(sub="sub", scopes=["projects"]):

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
        project_slug = project["slug"]
        dataset_id = project["dataset"][0]["id"]

        # ADD RESOURCE TO DATASET

        res = await client.post(
            f"{PREFIX}/project/{project_slug}/{dataset_id}",
            json={
                "path": testdata_path.as_posix(),
                "type": "png",
            },
        )
        assert res.status_code == 201

        # ADD GLOBAL TASKS

        task_id_list = []
        res = await client.post(
            f"{PREFIX}/task/",
            json=dict(
                name="Create ZARR Structure",
                input_type="png",
                output_type="zarr",
            ),
        )
        assert res.status_code == 201
        task_id_list.append(res.json()["id"])
        res = await client.post(
            f"{PREFIX}/task/",
            json=dict(
                name="Yokogawa to ZARR", input_type="zarr", output_type="zarr"
            ),
        )
        assert res.status_code == 201
        task_id_list.append(res.json()["id"])

        # ADD GLOBAL WORKFLOW

        res = await client.post(
            f"{PREFIX}/workflow/",
            json=dict(
                name="test workflow",
            ),
        )
        assert res.status_code == 201
        workflow_id = res.json()["id"]

        # ADD TASK TO GLOBAL WORKFLOW

        for task_id in task_id_list:
            res = await client.post(
                f"{PREFIX}/workflow/{workflow_id}",
                json=dict(
                    task_id=task_id,
                ),
            )
            assert res.status_code == 201

        # EXECUTE WORKFLOW

        res = await client.post(
            f"{PREFIX}/project/apply/{project_slug}"
            f"/{dataset_id}/{workflow_id}"
        )
        assert res.status_code == 201
