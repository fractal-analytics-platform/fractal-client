from pathlib import Path

from devtools import debug

from fractal_server.app.api.v1.task import _background_collect_pip
from fractal_server.app.api.v1.task import _TaskCollectPip
from fractal_server.app.api.v1.task import create_package_dir_pip
from fractal_server.app.api.v1.task import TaskCollectStatus
from fractal_server.app.models import State
from fractal_server.config import get_settings
from fractal_server.syringe import Inject
from fractal_server.tasks.collection import get_collection_path
from fractal_server.tasks.collection import get_log_path


async def test_task_get_list(db, client, task_factory, MockCurrentUser):
    t0 = await task_factory(name="task0")
    t1 = await task_factory(name="task1")
    t2 = await task_factory(index=2, subtask_list=[t0, t1])

    async with MockCurrentUser(persist=True):
        res = await client.get("api/v1/task/")
        data = res.json()
        assert res.status_code == 200
        debug(data)
        assert len(data) == 3
        assert data[2]["id"] == t2.id


async def test_background_collection(db, dummy_task_package):
    """
    GIVEN a package and its installation environment
    WHEN the background collection is called on it
    THEN the tasks are collected and the state is updated to db accordingly
    """
    task_pkg = _TaskCollectPip(package=dummy_task_package.as_posix())
    venv_path = create_package_dir_pip(
        task_pkg=task_pkg, user="test_bg_collection"
    )
    collection_status = TaskCollectStatus(
        status="pending", venv_path=venv_path, package=task_pkg.package
    )
    # replacing with path because of non-serializable Path
    collection_status_dict = collection_status.sanitised_dict()

    state = State(data=collection_status_dict)
    db.add(state)
    await db.commit()
    await db.refresh(state)
    debug(state)
    tasks = await _background_collect_pip(
        state=state, venv_path=venv_path, task_pkg=task_pkg, db=db
    )
    debug(tasks)
    assert tasks
    out_state = await db.get(State, state.id)
    debug(out_state)
    assert out_state.data["status"] == "OK"


async def test_collection_api(client, dummy_task_package, MockCurrentUser):
    """
    GIVEN a package in a format that `pip` understands
    WHEN the api to collect tasks from that package is called
    THEN
        * a dedicated directory is created and returned
        * in the background, an environment is created, the package is
          installed and the task collected
        * it is possible to GET the collection with the path to the folder to
          check the status of the background process
    """
    PREFIX = "/api/v1/task"

    task_collection = dict(package=dummy_task_package.as_posix())

    async with MockCurrentUser(persist=True):
        # NOTE: collecting private tasks so that they are assigned to user and
        # written in a non-default folder. Bypass for non stateless
        # FRACTAL_ROOT in test suite.
        res = await client.post(
            f"{PREFIX}/collect/pip/?public=false", json=task_collection
        )
        debug(res.json())
        assert res.status_code == 201
        assert res.json()["data"]["status"] == "pending"

        state = res.json()
        data = state["data"]
        assert "fractal_tasks_dummy" in data["venv_path"]
        venv_path = Path(data["venv_path"])

        res = await client.get(f"{PREFIX}/collect/{state['id']}")
        debug(res.json())
        assert res.status_code == 200
        state = res.json()
        data = state["data"]

        assert data["status"] == "OK"
        task_list = data["task_list"]
        assert data["log"] is None

        task_names = (t["name"] for t in task_list)
        assert len(task_list) == 2
        assert "dummy" in task_names
        assert "dummy parallel" in task_names

        # using verbose option
        res = await client.get(f"{PREFIX}/collect/{state['id']}?verbose=true")
        debug(res.json())
        state = res.json()
        data = state["data"]
        assert res.status_code == 200
        assert data["log"] is not None

        settings = Inject(get_settings)
        full_path = settings.FRACTAL_ROOT / venv_path
        assert get_collection_path(full_path).exists()
        assert get_log_path(full_path).exists()


async def test_collection_api_invalid_manifest(
    client, dummy_task_package_invalid_manifest, MockCurrentUser
):
    """
    GIVEN a package in a format that `pip` understands, with invalid manifest
    WHEN the api to collect tasks from that package is called
    THEN it returns 422 (Unprocessable Entity)
    """
    PREFIX = "/api/v1/task"

    task_collection = dict(
        package=dummy_task_package_invalid_manifest.as_posix()
    )

    async with MockCurrentUser(persist=True):
        # NOTE: collecting private tasks so that they are assigned to user and
        # written in a non-default folder. Bypass for non stateless
        # FRACTAL_ROOT in test suite.
        res = await client.post(
            f"{PREFIX}/collect/pip/?public=false", json=task_collection
        )
        debug(res.json())
        assert res.status_code == 422
