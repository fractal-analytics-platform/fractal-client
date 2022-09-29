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
from typing import Any
from typing import Dict
from typing import List

from devtools import debug


PREFIX = "/api/v1"


def task_id_by_name(name: str, task_list: List[Dict[str, Any]]) -> int:
    for task in task_list:
        if task["name"] == name:
            return task["id"]
    raise ValueError("No task `{name}`")


async def test_full_workflow(
    client,
    MockCurrentUser,
    testdata_path,
    collect_tasks,
    tmp_path,
):
    async with MockCurrentUser(persist=True):

        # CREATE PROJECT
        res = await client.post(
            f"{PREFIX}/project/",
            json=dict(
                name="test project",
                project_dir=tmp_path.as_posix(),
            ),
        )
        assert res.status_code == 201
        project = res.json()
        project_id = project["id"]
        input_dataset_id = project["dataset_list"][0]["id"]

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
                type="zarr",
            ),
        )
        debug(res.json())
        assert res.status_code == 201
        output_dataset = res.json()
        output_dataset_id = output_dataset["id"]

        res = await client.post(
            f"{PREFIX}/project/{project_id}/{output_dataset['id']}",
            json=dict(path=tmp_path.as_posix(), glob_pattern="*.zarr"),
        )
        out_resource = res.json()
        debug(out_resource)
        assert res.status_code == 201

        # CHECK WHERE WE ARE AT
        res = await client.get(f"{PREFIX}/project/{project_id}")
        debug(res.json())

        # CREATE WORKFLOW

        res = await client.post(
            f"{PREFIX}/task/",
            json=dict(
                name="my workflow",
                resource_type="workflow",
                input_type="image",
                output_type="zarr",
            ),
        )
        wf = res.json()
        workflow_id = wf["id"]
        debug(wf)
        assert res.status_code == 201

        res = await client.get(f"{PREFIX}/task/")
        assert res.status_code == 200
        task_list = res.json()

        task_id_create_zarr = task_id_by_name(
            name="Create OME-ZARR structure", task_list=task_list
        )
        task_id_yokogawa = task_id_by_name(
            name="Yokogawa to Zarr", task_list=task_list
        )

        # add subtasks
        res = await client.post(
            f"{PREFIX}/task/{workflow_id}/subtask/",
            json=dict(
                subtask_id=task_id_create_zarr,
                args=dict(channel_parameters={"A01_C01": {}}),
            ),
        )
        assert res.status_code == 201
        res = await client.post(
            f"{PREFIX}/task/{workflow_id}/subtask/",
            json=dict(
                subtask_id=task_id_yokogawa,
                args=dict(parallelization_level="well"),
            ),
        )
        assert res.status_code == 201
        debug(res.json())

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
        res = await client.get(f"{PREFIX}/project/{project_id}")
        debug(res.json())

        out_ds = [
            ds
            for ds in res.json()["dataset_list"]
            if ds["name"] == "output dataset"
        ][0]

        output_path = out_ds["resource_list"][0]["path"]
        component = out_ds["meta"]["well"][0]

        debug(output_path)
        zarrurl = output_path + "/" + component + "0"
        debug(zarrurl)
        try:
            import dask.array as da

            data_czyx = da.from_zarr(zarrurl)
            assert data_czyx.shape == (1, 2, 2160, 2 * 2560)
            assert data_czyx[0, 0, 0, 0].compute() == 0
        except ImportError:
            pass


async def test_full_workflow_repeated_tasks(
    app,
    client,
    MockCurrentUser,
    collect_tasks,
    tmp_path,
):

    num_subtasks = 3

    async with MockCurrentUser(persist=True):

        # CREATE PROJECT
        res = await client.post(
            f"{PREFIX}/project/",
            json=dict(
                name="test project",
                project_dir=tmp_path.as_posix(),
            ),
        )
        assert res.status_code == 201
        project = res.json()
        project_id = project["id"]

        # ADD A RESOURCE TO THE INPUT/OUTPUT DATASET
        input_dataset_id = project["dataset_list"][0]["id"]
        res = await client.post(
            f"{PREFIX}/project/{project_id}/{input_dataset_id}",
            json={
                "path": tmp_path.as_posix(),
                "glob_pattern": "*.json",
            },
        )
        debug(res.json())
        assert res.status_code == 201

        # CHECK WHERE WE ARE AT
        res = await client.get(f"{PREFIX}/project/{project_id}")
        debug(res.json())

        # CREATE WORKFLOW
        res = await client.post(
            f"{PREFIX}/task/",
            json=dict(
                name="my workflow",
                resource_type="workflow",
                input_type="Any",
                output_type="None",
            ),
        )
        wf = res.json()
        workflow_id = wf["id"]
        debug(wf)
        assert res.status_code == 201

        # Extract ID of task "dummy"
        res = await client.get(f"{PREFIX}/task/")
        assert res.status_code == 200
        task_list = res.json()
        task_id = task_id_by_name(name="dummy", task_list=task_list)

        # add subtasks
        for ind_task in range(num_subtasks):
            res = await client.post(
                f"{PREFIX}/task/{workflow_id}/subtask/",
                json=dict(
                    subtask_id=task_id,
                    args=dict(iteration=ind_task),
                ),
            )
            assert res.status_code == 201

        # EXECUTE WORKFLOW

        payload = dict(
            project_id=project_id,
            input_dataset_id=input_dataset_id,
            output_dataset_id=input_dataset_id,
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
