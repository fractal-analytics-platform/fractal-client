import pytest
from fractal_server.app.security import FRACTAL_DEFAULT_GROUP_NAME

EMAIL_USER = "test@testmail.com"
PWD_USER = "testpassword"


def test_group_commands_auth(register_user, invoke, caplog):
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


def test_group_commands(user_factory, invoke_as_superuser):

    # get default group id and superuser id
    res = invoke_as_superuser("group list --user-ids")
    assert res.retcode == 0
    assert len(res.data) == 1  # Only one group (default)
    assert res.data[0]["name"] == FRACTAL_DEFAULT_GROUP_NAME
    assert len(res.data[0]["user_ids"]) == 1  # Only one user (superuser)

    default_group_id = res.data[0]["id"]
    superuser_id = res.data[0]["user_ids"][0]

    # create 3 standard users (by default in default group)
    user1 = user_factory(email="user1@fractal.xy", password="psw1")
    user1_id = user1["id"]
    assert user1["group_ids"] == [default_group_id]
    user2 = user_factory(email="user2@fractal.xy", password="psw2")
    user2_id = user2["id"]
    assert user2["group_ids"] == [default_group_id]
    user3 = user_factory(email="user3@fractal.xy", password="psw3")
    user3_id = user3["id"]
    assert user3["group_ids"] == [default_group_id]

    res = invoke_as_superuser("group list --user-ids")
    assert len(res.data) == 1
    assert res.data[0]["user_ids"] == [
        superuser_id,
        user1_id,
        user2_id,
        user3_id,
    ]

    # Create 2 new empty groups (`group new`)

    with pytest.raises(SystemExit):
        # missing 'name'
        invoke_as_superuser("group new")

    res = invoke_as_superuser("group new foo")
    assert res.retcode == 0
    assert res.data["name"] == "foo"
    assert res.data["user_ids"] == []
    group1_id = res.data["id"]

    res = invoke_as_superuser("group new bar")
    group2_id = res.data["id"]

    # Add users to groups (`group update`)

    with pytest.raises(SystemExit):
        # missing 'group_id' and 'new_user_ids'
        invoke_as_superuser("group update")
    with pytest.raises(SystemExit):
        # missing 'new_user_ids'
        invoke_as_superuser(f"group update {default_group_id}")
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

    # add `user3` and `user2` to `group2`
    res = invoke_as_superuser(
        f"group update {group2_id} --new-user-ids {user3_id} {user2_id}"
    )
    assert res.retcode == 0
    assert res.data["id"] == group2_id
    assert res.data["user_ids"] == [user2_id, user3_id]  # order
    # add also `superuser` to `group2`
    res = invoke_as_superuser(
        f"group update {group2_id} --new-user-ids {superuser_id}"
    )
    assert res.data["user_ids"] == [superuser_id, user2_id, user3_id]  # order

    # Check groups are updated

    res = invoke_as_superuser("group list --user-ids")
    assert len(res.data) == 3
    assert res.data[0]["id"] == default_group_id
    assert res.data[0]["user_ids"] == [
        superuser_id,
        user1_id,
        user2_id,
        user3_id,
    ]
    assert res.data[1]["id"] == group1_id
    assert res.data[1]["user_ids"] == [user1_id, user2_id]
    assert res.data[2]["id"] == group2_id
    # sorted in order of insertion
    assert res.data[2]["user_ids"] == [user3_id, user2_id, superuser_id]

    # Test `group get` command

    with pytest.raises(SystemExit):
        # missing 'group_id'
        invoke_as_superuser("group get")

    res = invoke_as_superuser(f"group get {default_group_id}")
    assert res.retcode == 0
    assert res.data["name"] == FRACTAL_DEFAULT_GROUP_NAME
    assert len(res.data["user_ids"]) == 4

    # Test `list` without `--user-ids`

    res = invoke_as_superuser("group list")
    for group in res.data:
        assert group["user_ids"] is None

    # Test `--batch`

    res = invoke_as_superuser("--batch group list")
    assert res.data == f"{default_group_id} {group1_id} {group2_id}"

    res = invoke_as_superuser(f"--batch group get {group1_id}")
    assert res.data == group1_id

    res = invoke_as_superuser("--batch group new xyz")
    assert isinstance(res.data, int)

    res = invoke_as_superuser(
        f"--batch group update {group2_id} --new-user-ids {user1_id}"
    )
    assert res.data == group2_id
