from pathlib import Path

import pytest
from devtools import debug  # noqa

from fractal_server.app.models import Task
from fractal_server.app.models import Workflow
from fractal_server.app.models import WorkflowTask
from fractal_server.app.runner._process import process_workflow
from fractal_server.tasks import dummy as dummy_module

PROJECT_ID = 1
WORKFLOW_ID_1 = 11
WORKFLOW_ID_2 = 12
INDEX1 = 1
INDEX2 = 2
MESSAGE1 = "This is message 1"
MESSAGE2 = "This is message 2"
DUMMY_TASK_NAME = "dummy task"
DEFAULT_ARGS = dict(message="Nothing to say", raise_error=False)

# Define Task dummy_task
dummy_task = Task(
    project_id=PROJECT_ID,
    default_args=DEFAULT_ARGS,
    command=f"python {dummy_module.__file__}",
    name=DUMMY_TASK_NAME,
)

# Define Workflow made of a single WorkflowTask
workflow_1 = Workflow(
    id=WORKFLOW_ID_1,
    task_list=[
        WorkflowTask(
            id=101,
            workflow_id=WORKFLOW_ID_1,
            task_id=dummy_task.id,
            task=dummy_task,
            args=dict(message=MESSAGE1, index=INDEX1),
            order=0,
        ),
    ],
)

# Define Workflow made of two WorkflowTask's
workflow_2 = Workflow(
    id=WORKFLOW_ID_2,
    task_list=[
        WorkflowTask(
            id=101,
            workflow_id=WORKFLOW_ID_2,
            task_id=dummy_task.id,
            task=dummy_task,
            args=dict(message=MESSAGE1, index=INDEX1),
            order=0,
        ),
        WorkflowTask(
            id=102,
            workflow_id=WORKFLOW_ID_2,
            task_id=dummy_task.id,
            task=dummy_task,
            args=dict(message=MESSAGE2, index=INDEX2),
            order=1,
        ),
    ],
)


@pytest.mark.parametrize("workflow", [workflow_1, workflow_2])
async def test_process_workflow(workflow: Workflow, tmp_path: Path):

    # Make all folders workflow-specific:
    WORKFLOW_PATH = tmp_path / f"workflow_dir_{workflow.id}"
    INPUT_PATH = tmp_path / f"input_{workflow.id}"
    # FIXME: what should go in OUTPUT_PATH? a folder or a file?
    OUTPUT_PATH = tmp_path / f"output_{workflow.id}"

    INPUT_PATH.mkdir(exist_ok=True, parents=True)
    OUTPUT_PATH.mkdir(exist_ok=True, parents=True)
    WORKFLOW_PATH.mkdir(exist_ok=True, parents=True)
    LOGGER_NAME = "test logger"

    debug(workflow)
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
    assert output_task_pars.metadata["history"] == [DUMMY_TASK_NAME] * len(
        workflow.task_list
    )
