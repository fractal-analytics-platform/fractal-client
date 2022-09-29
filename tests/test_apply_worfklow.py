import json
from copy import deepcopy
from pathlib import Path

import pytest
from devtools import debug
from sqlmodel import select

from fractal_server.app.models import Dataset
from fractal_server.app.models import Subtask
from fractal_server.app.models import Task
from fractal_server.app.models import TaskRead


LEN_NONTRIVIAL_WORKFLOW = 3


@pytest.fixture
def nontrivial_workflow():
    workflow = Task(
        name="outer workflow",
        resource_type="workflow",
        subtask_list=[
            Subtask(
                subtask=Task(
                    name="inner workflow",
                    resource_type="workflow",
                    subtask_list=[
                        Subtask(
                            args={"message": "dummy0"},
                            subtask=Task(
                                name="dummy0",
                                module="fractal_tasks_core.dummy:dummy",
                                default_args=dict(
                                    message="dummy0", executor="cpu-low"
                                ),
                            ),
                        ),
                        Subtask(
                            args={"message": "dummy1"},
                            subtask=Task(
                                name="dummy1",
                                module="fractal_tasks_core.dummy:dummy",
                                default_args=dict(
                                    message="dummy1", executor="cpu-low"
                                ),
                            ),
                        ),
                    ],
                ),
            ),
            Subtask(
                subtask=Task(
                    name="dummy2",
                    module="fractal_tasks_core.dummy:dummy",
                    default_args=dict(message="dummy2", executor="cpu-low"),
                )
            ),
        ],
    )
    return workflow


N_INDICES = 3
DUMMY_MESSAGE = "dummy task message"
DUMMY_SUBTASK_MESSAGE = "dummy subtask message"


dummy_task = Task(
    name="dummy",
    resource_type="core task",
    module="fractal_tasks_core.dummy:dummy",
    default_args={"message": DUMMY_MESSAGE, "executor": "cpu-low"},
)
dummy_subtask = Subtask(
    subtask=dummy_task, args={"message": DUMMY_SUBTASK_MESSAGE}
)
dummy_subtask_parallel = Subtask(
    subtask=dummy_task,
    args={"message": DUMMY_SUBTASK_MESSAGE, "parallelization_level": "index"},
)


@pytest.mark.parametrize(
    ("task", "message", "nfiles"),
    [
        (dummy_task, DUMMY_MESSAGE, 1),
        (dummy_subtask, DUMMY_SUBTASK_MESSAGE, 1),
        (dummy_subtask_parallel, DUMMY_SUBTASK_MESSAGE, N_INDICES),
    ],
)
def test_atomic_task_factory(task, message, nfiles, tmp_path, patch_settings):
    """
    GIVEN
        * a task or subtask
        * input_path iterable
        * output_path
    WHEN passed to the task factory
    THEN
        * the relative PARSL workflow is correctly
        * it can run
        * the output is as expected
    """

    from fractal_server.app.runner import _atomic_task_factory
    from fractal_server.app.runner.runner_utils import load_parsl_config

    workflow_id = 0

    load_parsl_config(enable_monitoring=False, workflow_id=workflow_id)

    input_path_str = "/input/path"
    output_path = tmp_path
    metadata = {"index": list(range(N_INDICES))}

    parsl_app_future = _atomic_task_factory(
        task=task,
        input_paths=[Path(input_path_str)],
        output_path=output_path,
        metadata=metadata,
        workflow_id=workflow_id,
    )

    debug(parsl_app_future)
    metadata = parsl_app_future.result()
    debug(metadata)
    assert metadata

    assert sum(1 for item in output_path.glob("*.json")) == nfiles

    for r in output_path.glob("*.json"):
        with open(r, "r") as output_file:
            data = json.load(output_file)
            debug(data)
            assert len(data) == 1
            assert data[0]["message"] == message


def test_preprocess_workflow(nontrivial_workflow):
    """
    GIVEN a workflow with nested tasks
    WHEN the workflow is preprocessed
    THEN
        * the workflow is correctly unwrapped into a list
        * the order of the tasks is correctly preserved
    """
    workflow = nontrivial_workflow

    debug(workflow.preprocess())
    preprocessed_workflow = workflow.preprocess()
    for i, preprocessed_task in enumerate(preprocessed_workflow):
        assert str(i) in preprocessed_task.name

    assert i + 1 == LEN_NONTRIVIAL_WORKFLOW


def test_process_workflow(tmp_path, nontrivial_workflow, patch_settings):
    """
    GIVEN a nontrivial workflow
    WHEN the workflow is processed
    THEN
        * a single PARSL python_app which will execute the workflow is produced
        * it is executable
        * the output is the one expected from the workflow
    """

    from fractal_server.app.runner import _process_workflow

    app = _process_workflow(
        task=nontrivial_workflow,
        input_paths=[tmp_path / "0.json"],
        output_path=tmp_path / "0.json",
        metadata={},
    )
    debug(app)
    app.result()

    print(list(tmp_path.glob("*.json")))
    for f in tmp_path.glob("*.json"):
        with open(f, "r") as output_file:
            data = json.load(output_file)
            debug(data)
    assert len(data) == LEN_NONTRIVIAL_WORKFLOW
    assert data[0]["message"] == "dummy0"
    assert data[1]["message"] == "dummy1"
    assert data[2]["message"] == "dummy2"


def test_process_workflow_with_wrong_executor(tmp_path, patch_settings):
    """
    GIVEN a trivial workflow, with an invalid executor
    WHEN the workflow is processed
    THEN ValueError
    """

    from fractal_server.app.runner import _process_workflow

    dummy_task = Task(
        name="dummy",
        resource_type="core task",
        module="fake_module.dummy:dummy",
        default_args={"executor": "WRONG EXECUTOR"},
    )

    app = _process_workflow(
        task=dummy_task,
        input_paths=[tmp_path / "0.json"],
        output_path=tmp_path / "0.json",
        metadata={},
    )
    debug(app)

    with pytest.raises(ValueError):
        app.result()


async def test_apply_workflow(
    db,
    client,
    collect_tasks,
    MockCurrentUser,
    project_factory,
    dataset_factory,
    resource_factory,
    task_factory,
    tmp_path,
    patch_settings,
):
    """
    GIVEN
        * an input dataset and relative resource(s)
        * an output dataset and relative resource
        * a non-trivial workflow
    WHEN one applys the workflow to the input dataset
    THEN
        * the workflow is executed correctly
        * the output is correctly written in the output resource
    """

    from fractal_server.app.runner import submit_workflow

    # CREATE RESOURCES
    async with MockCurrentUser(persist=True) as user:
        prj = await project_factory(user)
        ds = await dataset_factory(prj, type="image")
        out_ds = await dataset_factory(prj, type="image", name="out_ds")

        await resource_factory(ds)
        output_path = (tmp_path / "0.json").as_posix()
        await resource_factory(out_ds, path=output_path, glob_pattern=None)

    # CREATE NONTRIVIAL WORKFLOW
    wf = await task_factory(
        name="worfklow",
        module=None,
        resource_type="workflow",
        input_type="image",
    )

    stm = select(Task).where(Task.name == "dummy")
    res = await db.execute(stm)
    dummy_task = res.scalar()

    MESSAGE = "test apply workflow"
    await wf.add_subtask(db, subtask=dummy_task, args=dict(message=MESSAGE))
    debug(TaskRead.from_orm(wf))

    # DONE CREATING WORKFLOW

    await submit_workflow(
        db=db, input_dataset=ds, output_dataset=out_ds, workflow=wf
    )
    with open(output_path, "r") as f:
        data = json.load(f)
        debug(data)
    assert len(data) == 1
    assert data[0]["message"] == MESSAGE

    await db.refresh(out_ds)
    debug(out_ds)
    assert out_ds.meta


async def test_create_zarr(
    db,
    client,
    collect_tasks,
    MockCurrentUser,
    project_factory,
    dataset_factory,
    resource_factory,
    task_factory,
    tmp_path,
    patch_settings,
):
    """
    GIVEN
        * some test png images
        * create ome-zarr structure task
        * a project, dataset and resource that represent the images
    WHEN
        * the task is applied on the resource
    THEN
        * the ZARR structure is correctly created
    """

    from fractal_server.app.runner import submit_workflow

    # CREATE RESOURCES
    async with MockCurrentUser(persist=True) as user:
        prj = await project_factory(user)
        ds = await dataset_factory(prj, type="image")
        out_ds = await dataset_factory(prj, type="zarr", name="out_ds")

        await resource_factory(ds)
        output_path = (tmp_path / "*.zarr").as_posix()
        await resource_factory(out_ds, path=output_path, glob_pattern=None)

    # CREATE NONTRIVIAL WORKFLOW
    wf = await task_factory(
        name="worfklow",
        module=None,
        resource_type="workflow",
        input_type="image",
    )

    stm = select(Task).where(Task.name == "Create OME-ZARR structure")
    res = await db.execute(stm)
    create_ome_zarr_task = res.scalar()

    # Modify the task default args
    db_task = await db.get(Task, create_ome_zarr_task.id)
    current_default_args = deepcopy(db_task._arguments)
    current_default_args.update(dict(channel_parameters={"A01_C01": {}}))
    setattr(db_task, "default_args", current_default_args)
    await db.commit()
    await db.refresh(db_task)

    await wf.add_subtask(db, subtask=create_ome_zarr_task)
    debug(TaskRead.from_orm(wf))

    # DONE CREATING WORKFLOW

    await submit_workflow(
        db=db, input_dataset=ds, output_dataset=out_ds, workflow=wf
    )
    zattrs = Path(output_path).parent / "myplate.zarr/.zattrs"
    with open(zattrs) as f:
        data = json.load(f)
        debug(data)
    assert len(data["plate"]["wells"]) == 1

    out_ds = await db.get(Dataset, out_ds.id)
    debug(out_ds)
    assert out_ds.meta
    # FIXME
    # The assertion above needs be specified to the metadata of the task


async def test_yokogawa(
    db,
    client,
    collect_tasks,
    MockCurrentUser,
    project_factory,
    dataset_factory,
    resource_factory,
    task_factory,
    tmp_path,
    patch_settings,
):
    """
    GIVEN
        * some test png images
        * create ome-zarr structure task
        * a project, dataset and resource that represent the images
    WHEN
        * the task is applied on the resource
    THEN
        * the ZARR structure is correctly created
    """

    from fractal_server.app.runner import submit_workflow

    # CREATE RESOURCES
    async with MockCurrentUser(persist=True) as user:
        prj = await project_factory(user)
        ds = await dataset_factory(prj, type="image")
        out_ds = await dataset_factory(prj, type="image", name="out_ds")

        await resource_factory(ds)
        output_path = (tmp_path / "*.zarr").as_posix()
        await resource_factory(out_ds, path=output_path, glob_pattern=None)

    # CREATE NONTRIVIAL WORKFLOW
    wf = await task_factory(
        name="worfklow",
        module=None,
        resource_type="workflow",
        input_type="image",
    )

    stm = select(Task).where(Task.name == "Create OME-ZARR structure")
    res = await db.execute(stm)
    create_ome_zarr_task = res.scalar()

    # Modify the task default args
    db_task = await db.get(Task, create_ome_zarr_task.id)
    current_default_args = deepcopy(db_task._arguments)
    current_default_args.update(dict(channel_parameters={"A01_C01": {}}))
    setattr(db_task, "default_args", current_default_args)
    await db.commit()
    await db.refresh(db_task)

    await wf.add_subtask(
        db,
        subtask=create_ome_zarr_task,
    )
    debug(TaskRead.from_orm(wf))

    stm = select(Task).where(Task.name == "Yokogawa to Zarr")
    res = await db.execute(stm)
    yokogawa = res.scalar()
    await wf.add_subtask(
        db, subtask=yokogawa, args=dict(parallelization_level="well")
    )
    debug(TaskRead.from_orm(wf))

    # DONE CREATING WORKFLOW

    await submit_workflow(
        db=db, input_dataset=ds, output_dataset=out_ds, workflow=wf
    )
    out_ds = await db.get(Dataset, out_ds.id)
    debug(out_ds)

    assert (
        out_ds.meta["history"][-1]
        == "Yokogawa to Zarr: ['myplate.zarr/B/03/0/']"
    )
    debug(out_ds.resource_list[0])
    zarrurl = (
        Path(out_ds.resource_list[0].path).parent
        / out_ds.meta["well"][0]
        / "0"
    ).as_posix()
    debug(zarrurl)

    try:
        import dask.array as da

        data_czyx = da.from_zarr(zarrurl)
        assert data_czyx.shape == (1, 2, 2160, 2 * 2560)
        assert data_czyx[0, 0, 0, 0].compute() == 0
    except ImportError:
        pass
