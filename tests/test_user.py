import json

import pytest
from devtools import debug
from fractal_server.app.security import FRACTAL_DEFAULT_GROUP_NAME

PWD_USER = "1234"


def test_register_as_user(invoke, caplog):
    with pytest.raises(SystemExit):
        invoke("user register aaa bbb")
    debug(caplog.text)
    assert "403" in caplog.text


@pytest.mark.parametrize("is_superuser", [True, False])
def test_register_as_superuser(
    invoke_as_superuser, is_superuser: bool, new_name
):
    EMAIL_USER = f"{new_name()}@example.org"
    if is_superuser:
        res = invoke_as_superuser(
            f"user register {EMAIL_USER} {PWD_USER} --superuser"
        )
        debug(res.data)
        assert res.retcode == 0
        assert res.data["is_superuser"]
    else:
        res = invoke_as_superuser(
            f"user register {EMAIL_USER} {PWD_USER} "
            "--slurm-user SOMETHING --username X"
        )
        debug(res.data)
        assert res.retcode == 0
        assert not res.data["is_superuser"]
    assert res.data["email"] == EMAIL_USER

    # Test that new user is verified (note: for the moment we don't expose the
    # possibility of registering a non-verified user)
    assert res.data["is_verified"]


def test_register_with_ssh_settings(invoke_as_superuser, new_name, tmp_path):

    EMAIL_USER = f"{new_name()}@example.org"

    with pytest.raises(SystemExit, match="File does not exist"):
        invoke_as_superuser(
            f"user register {EMAIL_USER} {PWD_USER} --ssh-settings-json xy.z"
        )

    invalid_json = tmp_path / "not-a-json.foo"
    with invalid_json.open("w") as f:
        f.write("hello world")
    with pytest.raises(SystemExit, match="not a valid JSON"):
        invoke_as_superuser(
            f"user register {EMAIL_USER} {PWD_USER} "
            f"--ssh-settings-json {invalid_json}"
        )

    invalid_key_json = tmp_path / "invalid-key.json"
    with invalid_key_json.open("w") as f:
        json.dump(dict(invalid="invalid"), f)
    with pytest.raises(SystemExit, match="Invalid key"):
        invoke_as_superuser(
            f"user register {EMAIL_USER} {PWD_USER} "
            f"--ssh-settings-json {invalid_key_json}"
        )

    valid_json = tmp_path / "ssh-config.json"
    with valid_json.open("w") as f:
        json.dump(
            dict(
                ssh_host="SSH_HOST",
                ssh_private_key_path="/SSH_PRIVATE_KEY_PATH",
            ),
            f,
        )
    PROJECT_DIR = "/somewhere/"
    res = invoke_as_superuser(
        f"user register {EMAIL_USER} {PWD_USER} "
        f"--project-dir {PROJECT_DIR} "
        f"--ssh-settings-json {valid_json}"
    )
    assert res.retcode == 0
    assert res.data["settings"]["project_dir"] == PROJECT_DIR
    assert res.data["settings"]["ssh_host"] == "SSH_HOST"
    assert res.data["settings"]["ssh_private_key_path"] == (
        "/SSH_PRIVATE_KEY_PATH"
    )


def test_register_as_superuser_with_batch(invoke_as_superuser, new_name):
    EMAIL_USER = f"{new_name()}@example.org"
    # Register a user with the --batch flag
    res = invoke_as_superuser(f"--batch user register {EMAIL_USER} {PWD_USER}")
    user_id = res.data
    debug(user_id)
    assert res.retcode == 0
    # Check that the user exists
    res = invoke_as_superuser(f"user show {user_id}")
    debug(res.data)
    assert res.data["email"] == EMAIL_USER
    assert res.retcode == 0


def test_list_as_user(invoke, caplog):
    with pytest.raises(SystemExit):
        invoke("user list")
    debug(caplog.text)
    assert "403" in caplog.text


def test_list_as_superuser(invoke_as_superuser, superuser, tester):
    res = invoke_as_superuser("user list")
    debug(res.data)
    assert res.retcode == 0
    list_emails = [user["email"] for user in res.data]
    debug(list_emails)
    assert superuser["email"] in list_emails
    assert tester["email"] in list_emails


def test_show_as_user(invoke, invoke_as_superuser, caplog, new_name):
    EMAIL_USER = f"{new_name()}@example.org"
    # Register a new user
    res = invoke_as_superuser(f"user register {EMAIL_USER} {PWD_USER}")
    user_id = res.data["id"]
    # Call fractal user show
    with pytest.raises(SystemExit):
        invoke(f"user show {user_id}")
    debug(caplog.text)
    assert "403" in caplog.text


def test_show_as_superuser(invoke_as_superuser, new_name):
    EMAIL_USER = f"{new_name()}@example.org"
    # Register a new user
    res = invoke_as_superuser(f"user register {EMAIL_USER} {PWD_USER}")
    user_id = res.data["id"]
    # Call fractal user show
    invoke_as_superuser(f"user show {user_id}")
    debug(res.data)
    assert res.retcode == 0
    assert res.data["email"] == EMAIL_USER


def test_edit_as_user(invoke, invoke_as_superuser, caplog, new_name):
    EMAIL_USER = f"{new_name()}@example.org"
    # Register a new user
    res = invoke_as_superuser(f"user register {EMAIL_USER} {PWD_USER}")
    user_id = res.data["id"]
    # Call fractal user edit
    with pytest.raises(SystemExit):
        res = invoke(
            f"user edit {user_id} "
            "--new-email email@something.xy --make-verified"
        )
    debug(caplog.text)
    assert "403" in caplog.text


@pytest.mark.parametrize("new_is_superuser", [True, False])
@pytest.mark.parametrize("new_is_verified", [True, False])
@pytest.mark.parametrize("new_is_non_verified", [True, False])
def test_edit_as_superuser(
    invoke_as_superuser,
    new_is_superuser,
    new_is_verified,
    new_is_non_verified,
    new_name,
):
    EMAIL_USER = f"{new_name()}@example.org"
    # Register a new user
    res = invoke_as_superuser(f"user register {EMAIL_USER} {PWD_USER}")
    assert res.retcode == 0
    user_id = res.data["id"]
    # Call fractal user edit
    NEW_EMAIL = f"{new_name()}@example.org"
    NEW_SLURM_USER = "new_slurm"
    NEW_USERNAME = "new_username"
    cmd = (
        f"user edit {user_id} "
        f"--new-email {NEW_EMAIL} "
        f"--new-password SOMETHING "
        f"--new-slurm-user {NEW_SLURM_USER} "
        f"--new-username {NEW_USERNAME} "
    )
    if new_is_superuser:
        cmd = f"{cmd} --make-superuser"
    if new_is_verified:
        cmd = f"{cmd} --make-verified"
    if new_is_non_verified:
        cmd = f"{cmd} --remove-verified"

    if new_is_verified and new_is_non_verified:
        with pytest.raises(SystemExit):
            invoke_as_superuser(cmd)
    elif new_is_verified or new_is_non_verified:
        res = invoke_as_superuser(cmd)
        assert res.retcode == 0
        assert res.data["email"] == NEW_EMAIL
        assert res.data["username"] == NEW_USERNAME
        assert res.data["is_superuser"] == new_is_superuser
        assert (
            res.data["is_verified"]
            if new_is_verified
            else not res.data["is_verified"]
        )
        assert res.data["settings"]["slurm_user"] == NEW_SLURM_USER
    else:
        res = invoke_as_superuser(cmd)
        assert res.retcode == 1
        assert res.data == (
            "Cannot use `--new-email` without `--make-verified` or "
            "`--remove-verified`"
        )

    # If the user was made a superuser, check that we can go back to normal
    # user
    if new_is_superuser:
        cmd = f"user edit {user_id} --remove-superuser"
        debug(cmd)
        res = invoke_as_superuser(cmd)
        debug(res.data)
        assert res.retcode == 0
        assert not res.data["is_superuser"]

    # If the user was made verified, check that we can go back to normal
    # user
    if new_is_non_verified:
        cmd = f"user edit {user_id} --make-verified"
        debug(cmd)
        res = invoke_as_superuser(cmd)
        debug(res.data)
        assert res.retcode == 0
        assert res.data["is_verified"]


def test_edit_user_settings(invoke_as_superuser, tmp_path, new_name):
    EMAIL_USER = f"{new_name()}@example.org"

    EMPTY_USER_SETTINGS = {
        "ssh_host": None,
        "ssh_username": None,
        "ssh_private_key_path": None,
        "ssh_tasks_dir": None,
        "ssh_jobs_dir": None,
        "slurm_user": None,
        "slurm_accounts": [],
        "project_dir": None,
    }
    SSH_HOST = "something.somewhere"
    SSH_PRIVATE_KEY_PATH = "/tmp/something.key"
    NEW_PROJECT_DIR = "/somewhere/else/"
    NEW_USER_SETTINGS = {
        "ssh_host": SSH_HOST,
        "ssh_username": None,
        "ssh_private_key_path": SSH_PRIVATE_KEY_PATH,
        "ssh_tasks_dir": None,
        "ssh_jobs_dir": None,
        "slurm_user": None,
        "slurm_accounts": [],
        "project_dir": NEW_PROJECT_DIR,
    }

    # Register a new user
    res = invoke_as_superuser(f"user register {EMAIL_USER} {PWD_USER}")
    assert res.retcode == 0
    user_id = res.data["id"]

    # Check empty user settings
    res = invoke_as_superuser(f"user show {user_id}")
    assert res.retcode == 0
    user_settings = {
        key: value
        for key, value in res.data["settings"].items()
        if key != "id"
    }
    debug(user_settings)
    assert user_settings == EMPTY_USER_SETTINGS

    # Call fractal user edit
    ssh_settings_file = tmp_path / "ssh.json"
    with ssh_settings_file.open("w") as f:
        json.dump(
            {
                "ssh_host": SSH_HOST,
                "ssh_private_key_path": SSH_PRIVATE_KEY_PATH,
            },
            f,
        )
    cmd = (
        f"user edit {user_id} "
        f"--new-project-dir {NEW_PROJECT_DIR} "
        f"--new-ssh-settings-json {ssh_settings_file.as_posix()}"
    )
    res = invoke_as_superuser(cmd)
    assert res.retcode == 0
    debug(res.data)

    # Check edited user settings
    res = invoke_as_superuser(f"user show {user_id}")
    assert res.retcode == 0
    user_settings = {
        key: value
        for key, value in res.data["settings"].items()
        if key != "id"
    }
    debug(user_settings)
    assert user_settings == NEW_USER_SETTINGS

    # Failure due to missing file
    ssh_settings_file = tmp_path / "invalid-ssh.json"
    cmd = (
        f"user edit {user_id} "
        f"--new-ssh-settings-json {ssh_settings_file.as_posix()}"
    )
    with pytest.raises(SystemExit, match="File does not exist."):
        res = invoke_as_superuser(cmd)

    # Failure due to file not being a valid JSON
    invalid_json = tmp_path / "invalid-json.foo"
    with invalid_json.open("w") as f:
        f.write("hello world")
    cmd = (
        f"user edit {user_id} "
        f"--new-ssh-settings-json {invalid_json.as_posix()}"
    )
    with pytest.raises(SystemExit, match="not a valid JSON"):
        res = invoke_as_superuser(cmd)

    # Failure due to invalid keys
    ssh_settings_file = tmp_path / "invalid-ssh.json"
    with ssh_settings_file.open("w") as f:
        json.dump(
            dict(invalid="invalid"),
            f,
        )
    cmd = (
        f"user edit {user_id} "
        f"--new-ssh-settings-json {ssh_settings_file.as_posix()}"
    )
    with pytest.raises(SystemExit, match="Invalid key"):
        res = invoke_as_superuser(cmd)


def test_edit_arguments(invoke_as_superuser):
    # Test that superuser flags are mutually exclusive
    with pytest.raises(SystemExit):
        cmd = "user edit SOME_USER_ID --make-superuser --remove-superuser"
        invoke_as_superuser(cmd)


def test_whoami_as_user(invoke, tester):
    res = invoke("user whoami")
    assert res.retcode == 0
    debug(res.data)
    assert res.data["email"] == tester["email"]
    assert not res.data["is_superuser"]
    user_id = res.data["id"]

    # Test user whoami with --batch flag
    res = invoke("--batch user whoami")
    debug(res.data)
    assert res.data == user_id
    assert res.retcode == 0


def test_whoami_as_superuser(invoke_as_superuser, superuser):
    res = invoke_as_superuser("user whoami")
    assert res.retcode == 0
    debug(res.data)
    assert res.data["email"] == superuser["email"]
    assert res.data["is_superuser"]


def test_user_set_groups(invoke_as_superuser, user_factory, new_name):
    # get default group
    res = invoke_as_superuser("group list --user-ids")
    default_group = next(
        group
        for group in res.data
        if group["name"] == FRACTAL_DEFAULT_GROUP_NAME
    )
    default_group_id = default_group["id"]
    # create 2 new users
    user1 = user_factory(email=f"{new_name()}@example.org", password="psw1")
    assert len(user1["group_ids_names"]) == 1
    user1_id = user1["id"]
    user2 = user_factory(email=f"{new_name()}@example.org", password="psw2")
    assert len(user2["group_ids_names"]) == 1
    user2_id = user2["id"]
    # create 2 new groups
    group1 = invoke_as_superuser(f"group new {new_name()}")
    assert len(group1.data["user_ids"]) == 0
    group1_id = group1.data["id"]
    group2 = invoke_as_superuser(f"group new {new_name()}")
    assert len(group2.data["user_ids"]) == 0
    group2_id = group2.data["id"]

    with pytest.raises(SystemExit):
        # no arguments
        invoke_as_superuser("user set-groups")
    with pytest.raises(SystemExit):
        # no group_ids list
        invoke_as_superuser(f"user set-groups {user1_id}")
    with pytest.raises(SystemExit):
        # group_ids must be a list of integers
        invoke_as_superuser(f"user set-groups {user1_id} {group1_id} foo")
    with pytest.raises(SystemExit):
        # there must always be the default group id in group_ids
        invoke_as_superuser(f"user set-groups {user1_id} {group1_id}")
    with pytest.raises(SystemExit):
        # repeated elements in group_ids are forbidden
        invoke_as_superuser(
            "user set-groups "
            f"{user1_id} {default_group_id} {group1_id} {group1_id}"
        )

    # Add user1 to group1
    res = invoke_as_superuser(
        f"user set-groups {user1_id} {group1_id} {default_group_id}"
    )
    assert len(res.data["group_ids_names"]) == 2
    group1 = invoke_as_superuser(f"group get {group1_id}")
    assert len(group1.data["user_ids"]) == 1

    # Add user2 to group1 and group2
    res = invoke_as_superuser(
        "user set-groups "
        f"{user2_id} {group2_id} {group1_id} {default_group_id}"
    )
    assert len(res.data["group_ids_names"]) == 3
    group1 = invoke_as_superuser(f"group get {group1_id}")
    assert len(group1.data["user_ids"]) == 2
    group2 = invoke_as_superuser(f"group get {group2_id}")
    assert len(group2.data["user_ids"]) == 1

    # Add user1 to group2 and remove them from group1
    res = invoke_as_superuser(
        f"user set-groups {user1_id} {group2_id} {default_group_id}"
    )
    assert len(res.data["group_ids_names"]) == 2
    group1 = invoke_as_superuser(f"group get {group1_id}")
    assert len(group1.data["user_ids"]) == 1
    group2 = invoke_as_superuser(f"group get {group2_id}")
    assert len(group2.data["user_ids"]) == 2
