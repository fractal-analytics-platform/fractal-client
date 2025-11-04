import pytest
from fractal_client.cmd import dataset
from fractal_client.cmd import group
from fractal_client.cmd import job
from fractal_client.cmd import NoCommandError
from fractal_client.cmd import profile
from fractal_client.cmd import project
from fractal_client.cmd import resource
from fractal_client.cmd import task
from fractal_client.cmd import user
from fractal_client.cmd import workflow


def test_invalid_commands(invoke):
    for arg in ["", " INVALID"]:
        for command in [
            "",
            "dataset",
            "project",
            "job",
            "workflow",
            "task",
            "user",
            "group",
            "resource",
            "profile",
        ]:
            with pytest.raises(SystemExit):
                invoke(f"{command}{arg}")


def test_unit_invalid_subcommand():
    for _function in [
        project,
        dataset,
        task,
        workflow,
        job,
        user,
        group,
        resource,
        profile,
    ]:
        with pytest.raises(NoCommandError):
            _function(client=None, subcmd="invalid")
