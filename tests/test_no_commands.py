import pytest


def test_no_commands(invoke):
    for command in [
        "",
        "dataset",
        "project",
        "job",
        "workflow",
        "task",
        "user",
    ]:
        with pytest.raises(SystemExit):
            invoke(command)
