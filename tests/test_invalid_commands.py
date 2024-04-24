import pytest


def test_invalid_commands(invoke):
    for arg in ["", " INVALID"]:
        for command in [
            "",
            "dataset",
            "project",
            "job",
            "workflow",
            "task" "user",
        ]:
            with pytest.raises(SystemExit):
                invoke(f"{command}{arg}")
