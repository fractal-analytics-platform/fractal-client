import json
from pathlib import Path

import pytest
from devtools import debug

from fractal_client.cmd._aux_task_caching import TASKS_CACHE_FILENAME
from fractal_client.config import settings


COLLECTION_TIMEOUT = 15.0


def test_task_new(
    invoke,
    invoke_as_custom_user,
    tmp_path,
    new_name,
    user_factory,
):

    # create a new task with just positional required args
    args_path = str(tmp_path / "args.json")
    args = {"image_dir": "/asdasd"}
    with open(args_path, "w") as f:
        json.dump(args, f)

    meta_path = str(tmp_path / "meta.json")
    meta = {"a": "b"}
    with open(meta_path, "w") as f:
        json.dump(meta, f)

    TASK_NAME = new_name()
    res = invoke(
        f"task new {TASK_NAME} --command-parallel _command "
        f"--version _version --meta-parallel {meta_path} "
        f"--args-schema-parallel {args_path} "
        f"--args-schema-version 1.0.0 "
        "--private"
    )
    debug(res.data)
    assert res.retcode == 0
    assert res.data["name"] == TASK_NAME
    assert res.data["command_parallel"] == "_command"
    assert res.data["version"] == "_version"
    assert res.data["meta_parallel"] == meta
    assert res.data["args_schema_version"] == "1.0.0"
    first_task_id = int(res.data["id"])

    # Check that task is actually private
    new_user_credentials = dict(
        email=f"{new_name()}@example.org",
        password="1234",
    )
    user_factory(**new_user_credentials)
    with pytest.raises(SystemExit):
        res = invoke_as_custom_user(
            f"task show {first_task_id}",
            **new_user_credentials,
        )

    # create a new task with batch option
    TASK_NAME_2 = new_name()
    res = invoke(
        f"--batch task new {TASK_NAME_2} --command-parallel _command2"
    )
    res.show()
    assert res.retcode == 0
    assert res.data == str(first_task_id + 1)

    # create a new task with same name as before. Note that in check_response
    # we have sys.exit(1) when status code is not the expecte one
    with pytest.raises(SystemExit) as e:
        invoke(f"task new {TASK_NAME_2} --command-parallel _command2")
    assert e.value.code == 1

    # create a new task passing not existing file
    res = invoke(
        (
            f"task new {new_name()} --command-parallel _command "
            "--meta-parallel ./foo.pdf"
        )
    )
    assert res.retcode == 1

    metanp_path = str(tmp_path / "meta.json")
    metanp = {"a": "b"}
    with open(metanp_path, "w") as f:
        json.dump(metanp, f)
    res = invoke(
        (
            f"task new {new_name()} --command-non-parallel _command_np "
            f"--meta-non-parallel {metanp_path} "
            f"--args-schema-non-parallel {args_path} "
        )
    )
    assert res.data["args_schema_non_parallel"] == args


def test_task_edit(
    caplog,
    invoke,
    tmp_path,
    new_name,
):

    args_path = str(tmp_path / "args.json")
    args = {"image_dir": "/asdasd"}
    with open(args_path, "w") as f:
        json.dump(args, f)

    meta_path = str(tmp_path / "meta.json")
    meta = {"a": "b"}
    with open(meta_path, "w") as f:
        json.dump(meta, f)

    NAME = new_name()
    task = invoke(
        f"task new {NAME} --command-parallel _command "
        f"--version _version --meta-parallel {meta_path} "
        f"--args-schema-parallel {args_path} "
        f"--args-schema-version 1.0.0"
    )

    task.show()
    assert task.retcode == 0
    task_id = task.data["id"]

    # Test successful edit of string attributes
    NEW_COMMAND_PARALLEL = "run_parallel"
    res = invoke(
        (
            f"task edit --id {task_id} "
            f"--command-parallel {NEW_COMMAND_PARALLEL}"
        )
    )
    assert res.data["command_parallel"] == NEW_COMMAND_PARALLEL
    assert res.retcode == 0

    # Add non-parallel task and test command-non-parallel

    meta_path = str(tmp_path / "meta.json")
    meta = {"a": "b"}
    with open(meta_path, "w") as f:
        json.dump(meta, f)

    task_np = invoke(
        (
            f"task new {new_name()} --command-non-parallel _command_np "
            f"--version 1.0.1 --meta-non-parallel {meta_path}"
        )
    )

    NEW_COMMAND_NON_PARALLEL = "run_non_parallel"
    res = invoke(
        (
            f"task edit --id {task_np.data['id']} "
            f"--command-non-parallel {NEW_COMMAND_NON_PARALLEL}"
        )
    )
    assert res.data["command_non_parallel"] == NEW_COMMAND_NON_PARALLEL
    assert res.retcode == 0

    # Test fail with no task_id nor task_name
    with pytest.raises(SystemExit):
        res = invoke("task edit")
    # Test fail with both task_id and task_name
    with pytest.raises(SystemExit):
        res = invoke(f"task edit --id {task_id} --name {task.data['name']}")
    # Test fail with both task_id and task_version
    with pytest.raises(SystemExit):
        res = invoke(f"task edit --id {task_id} --version 1.2.3.4.5.6")
    assert caplog.records[-1].msg == (
        "Too many arguments: cannot provide both `id` and `version`."
    )
    # Test fail "name and wrong version"
    with pytest.raises(SystemExit):
        invoke("task delete --name INVALID_NAME --version INVALID_VERSION")

    input_types = {"input": True, "output": False}

    i_types_path = str(tmp_path / "itypes.json")
    with open(i_types_path, "w") as f:
        json.dump(input_types, f)

    output_types = {"input": False, "output": True}

    o_types_path = str(tmp_path / "otypes.json")
    with open(o_types_path, "w") as f:
        json.dump(output_types, f)

    # Test regular updates (both by id and name)
    res = invoke(f"task edit --id {task_id} --input-types {i_types_path}")
    assert res.data["input_types"] == input_types
    assert res.retcode == 0
    res = invoke(f"task edit --name {NAME} --output-types {o_types_path}")
    assert res.data["output_types"] == output_types
    assert res.retcode == 0

    # Test regular update by name, after deleting cache
    cache_dir = Path(settings.FRACTAL_CACHE_PATH)
    cache_file = cache_dir / TASKS_CACHE_FILENAME
    cache_file.unlink(missing_ok=True)

    res = invoke(f"task edit --name {NAME} --output-types {o_types_path}")
    assert res.data["output_types"] == output_types
    assert res.retcode == 0

    # Test failing invalid name
    fail_output_types = {"input": True, "output": False}

    f_o_types_path = str(tmp_path / "fotypes.json")
    with open(f_o_types_path, "w") as f:
        json.dump(fail_output_types, f)

    with pytest.raises(SystemExit):
        res = invoke(
            f"task edit --name INVALID_NAME --output-types {f_o_types_path}"
        )

    # Test regular update by name, after creating an invalid cache
    with cache_file.open("w") as f:
        json.dump([], f)

    new_output_types = {"input": False, "output": True}
    n_o_types_path = str(tmp_path / "notypes.json")
    with open(n_o_types_path, "w") as f:
        json.dump(new_output_types, f)

    res = invoke(f"task edit --name {NAME} --output-types {n_o_types_path}")
    assert res.data["output_types"] == new_output_types
    assert res.retcode == 0
