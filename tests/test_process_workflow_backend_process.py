from devtools import debug  # noqa

from fractal_server.app.models import Task
from fractal_server.app.models import Workflow
from fractal_server.app.models import WorkflowTask

# from fractal_server.app.runner._process import process_workflow


async def test_process_workflow(tmp_path):

    PROJECT_ID = 1
    WORKFLOW_ID = 123
    INDEX = 2
    MESSAGE = "This is a message"
    DEFAULT_ARGS = dict(message="Nothing to say", raise_error=False)

    task_1 = Task(project_id=PROJECT_ID, default_args=DEFAULT_ARGS)
    debug(task_1)
    workflow_task_1 = WorkflowTask(
        id=3,
        workflow_id=WORKFLOW_ID,
        task_id=task_1.id,
        task=task_1,
        args=dict(message=MESSAGE, index=INDEX),
    )
    debug(workflow_task_1)

    workflow = Workflow(id=WORKFLOW_ID, task_list=[workflow_task_1])
    debug(workflow)

    """
    LOGGER_NAME = "test logger"
    output_task_pars = await process_workflow(
        workflow=workflow,
        input_paths=[tmp_path / "input"],
        output_path=tmp_path / "output",
        input_metadata={},
        logger_name=LOGGER_NAME,
        workflow_dir=tmp_path/"workflow_dir",
        username=None,
    )

    debug(output_task_pars)
    """
