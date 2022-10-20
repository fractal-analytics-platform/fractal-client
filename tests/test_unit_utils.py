from datetime import timezone

from devtools import debug

from fractal_server.app.models.models_utils import get_timestamp


def test_timestamp():
    """
    GIVEN a function that provides a timestamp
    WHEN called
    THEN the timestamp is timezone aware and the timezone is set to UTC
    """
    ts = get_timestamp()
    debug(ts)
    assert ts
    assert ts.tzinfo is timezone.utc


def test_cascade_delete_workflow():
    """
    GIVEN a Workflow
    WHEN the Workflow is deleted
    THEN all the related WorkflowTask are deleted
    """
    pass


def test_cascade_delete_project():
    """
    GIVEN a Project
    WHEN the Project is deleted
    THEN all the related Workflows are deleted
    """
    pass
