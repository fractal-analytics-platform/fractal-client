import pytest
from fractal_server.app.security import FRACTAL_DEFAULT_GROUP_NAME


def test_group_commands_auth(invoke, caplog):
    """
    Assert 'group' commands are not accessible to standard users
    """

    def _assert_403(cmd):
        caplog.clear()
        with pytest.raises(SystemExit):
            invoke(cmd)
        assert "403" in caplog.text

    _assert_403(cmd="group list")
    _assert_403(cmd="group get 1")
    _assert_403(cmd="group new foo")
    _assert_403(cmd="group update 1 --new-user-ids 1")


def test_group_commands(
    user_factory, invoke_as_superuser, new_name, superuser
):

    # get default group id and superuser id
    res = invoke_as_superuser("group list --user-ids")
    assert res.retcode == 0
    initial_number_of_groups = len(res.data)
    default_group = next(
        group
        for group in res.data
        if group["name"] == FRACTAL_DEFAULT_GROUP_NAME
    )
    initial_number_of_users = len(default_group["user_ids"])

    default_group_id = default_group["id"]
    superuser_id = superuser["id"]

    # create 3 standard users (by default in default group)
    user1 = user_factory(email=f"{new_name()}@example.org", password="psw1")
    user1_id = user1["id"]
    assert user1["group_ids_names"] == [[default_group_id, "All"]]
    user2 = user_factory(email=f"{new_name()}@example.org", password="psw2")
    user2_id = user2["id"]
    assert user2["group_ids_names"] == [[default_group_id, "All"]]
    user3 = user_factory(email=f"{new_name()}@example.org", password="psw3")
    user3_id = user3["id"]
    assert user3["group_ids_names"] == [[default_group_id, "All"]]

    res = invoke_as_superuser("group list --user-ids")
    assert len(res.data) == initial_number_of_groups
    assert len(res.data[0]["user_ids"]) == initial_number_of_users + 3

    # Create 2 new empty groups (`group new`)

    with pytest.raises(SystemExit):
        # missing 'name'
        invoke_as_superuser("group new")

    NEW_NAME = new_name()
    res = invoke_as_superuser(f"group new {NEW_NAME} --viewer-paths /a /b")
    assert res.retcode == 0
    assert res.data["name"] == NEW_NAME
    assert res.data["user_ids"] == []
    group1_viewer_paths = res.data["viewer_paths"]
    assert group1_viewer_paths == ["/a", "/b"]
    group1_id = res.data["id"]

    res = invoke_as_superuser(f"group new {new_name()}")
    group2_id = res.data["id"]
    group2_viewer_paths = res.data["viewer_paths"]
    assert group2_viewer_paths == []

    # Add users to groups (`group update`)

    # empty update
    res = invoke_as_superuser(f"group update {default_group_id}")
    assert res.retcode == 0
    assert res.data["id"] == default_group_id

    with pytest.raises(SystemExit):
        # missing 'group_id' and 'new_user_ids'
        invoke_as_superuser("group update")
    with pytest.raises(SystemExit):
        # missing 'group_id'
        invoke_as_superuser(f"group update --new-user-ids {superuser_id}")
    with pytest.raises(SystemExit):
        # user already in group
        invoke_as_superuser(
            f"group update {default_group_id} --new-user-ids {superuser_id}"
        )
    with pytest.raises(SystemExit):
        # non existing user
        invoke_as_superuser(
            f"group update {default_group_id} --new-user-ids 9999"
        )

    # add `user1` and `user2` to `group1`
    res = invoke_as_superuser(
        f"group update {group1_id} --new-user-ids {user1_id} {user2_id}"
    )
    assert res.retcode == 0
    assert res.data["id"] == group1_id
    assert res.data["user_ids"] == [user1_id, user2_id]
    assert res.data["viewer_paths"] == group1_viewer_paths

    # add `user3` and `user2` to `group2`
    res = invoke_as_superuser(
        f"group update {group2_id} --new-user-ids {user3_id} {user2_id}"
    )
    assert res.retcode == 0
    assert res.data["id"] == group2_id
    assert set(res.data["user_ids"]) == set([user3_id, user2_id])
    # add also `superuser` to `group2`
    res = invoke_as_superuser(
        f"group update {group2_id} --new-user-ids {superuser_id}"
    )
    assert set(res.data["user_ids"]) == set([user3_id, user2_id, superuser_id])
    assert res.data["viewer_paths"] == group2_viewer_paths

    # Check groups are updated

    res = invoke_as_superuser("group list --user-ids")
    assert len(res.data) == initial_number_of_groups + 2
    assert res.data[0]["id"] == default_group_id
    assert len(res.data[0]["user_ids"]) == initial_number_of_users + 3
    assert res.data[-2]["id"] == group1_id
    assert set(res.data[1]["user_ids"]) == set([user1_id, user2_id])
    assert res.data[-1]["id"] == group2_id
    assert set(res.data[2]["user_ids"]) == set(
        [user3_id, user2_id, superuser_id]
    )

    # Test `group get` command

    with pytest.raises(SystemExit):
        # missing 'group_id'
        invoke_as_superuser("group get")

    res = invoke_as_superuser(f"group get {default_group_id}")
    assert res.retcode == 0
    assert res.data["name"] == FRACTAL_DEFAULT_GROUP_NAME
    assert len(res.data["user_ids"]) == initial_number_of_users + 3

    # Test `list` without `--user-ids`

    res = invoke_as_superuser("group list")
    for group in res.data:
        assert group["user_ids"] is None

    # Test `--batch`

    res = invoke_as_superuser("--batch group list")
    assert len(res.data.split(" ")) == initial_number_of_groups + 2

    res = invoke_as_superuser(f"--batch group new {new_name()}")
    assert isinstance(res.data, int)

    # Test update of viewer-paths

    res_pre_patch = invoke_as_superuser(f"group get {group1_id}")
    assert res_pre_patch.retcode == 0
    res_pre_patch.data.pop("viewer_paths")
    res_post_patch = invoke_as_superuser(
        f"group update {group1_id} --new-viewer-paths /a/b /c/d"
    )
    assert res_post_patch.retcode == 0
    viewer_paths_post_pach = res_post_patch.data.pop("viewer_paths")
    assert viewer_paths_post_pach == ["/a/b", "/c/d"]
    assert res_post_patch.data == res_pre_patch.data

    # Test `whoami --viewer-paths`

    invoke_as_superuser(
        f"group update {group1_id} --new-user-ids {superuser_id}"
    )
    assert "viewer_paths" not in superuser
    res = invoke_as_superuser("user whoami --viewer-paths")
    assert set(res.data.get("viewer_paths")) == {"/a/b", "/c/d"}
