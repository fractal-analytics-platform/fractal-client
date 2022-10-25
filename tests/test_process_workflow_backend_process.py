from pathlib import Path

from devtools import debug  # noqa

from fractal_server.app.models import Task
from fractal_server.app.models import Workflow
from fractal_server.app.models import WorkflowTask
from fractal_server.app.runner._process import process_workflow
from fractal_server.tasks import dummy as dummy_module


async def test_process_workflow(tmp_path):

    PROJECT_ID = 1
    WORKFLOW_ID = 123
    INDEX = 2
    MESSAGE = "This is a message"
    DEFAULT_ARGS = dict(message="Nothing to say", raise_error=False)
    WORKFLOW_PATH = tmp_path / "workflow_dir"
    INPUT_PATH = tmp_path / "input"
    OUTPUT_PATH = tmp_path / "output"

    INPUT_PATH.mkdir(exist_ok=True, parents=True)
    OUTPUT_PATH.mkdir(exist_ok=True, parents=True)
    Path(WORKFLOW_PATH).mkdir(exist_ok=True, parents=True)

    # Define task
    command = f"python {dummy_module.__file__}"
    task_1 = Task(
        project_id=PROJECT_ID, default_args=DEFAULT_ARGS, command=command
    )
    debug(task_1)

    # Define workflow task
    workflow_task_1 = WorkflowTask(
        id=3,
        workflow_id=WORKFLOW_ID,
        task_id=task_1.id,
        task=task_1,
        args=dict(message=MESSAGE, index=INDEX),
        order=0,
    )
    debug(workflow_task_1)

    # Define workflow
    workflow = Workflow(id=WORKFLOW_ID, task_list=[workflow_task_1])
    debug(workflow)

    LOGGER_NAME = "test logger"
    output_task_pars = await process_workflow(
        workflow=workflow,
        input_paths=[INPUT_PATH],
        output_path=OUTPUT_PATH,
        input_metadata={},
        logger_name=LOGGER_NAME,
        workflow_dir=WORKFLOW_PATH,
        username=None,
    )

    debug(output_task_pars)

    # FIXME validate output (either in output_task_pars or in the WORKFLOW_PATH
    # files)
