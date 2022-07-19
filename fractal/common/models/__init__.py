from typing import Optional

from pydantic import BaseModel


class ApplyWorkflow(BaseModel):
    project_id: int
    input_dataset_id: int
    output_dataset_id: Optional[int]
    workflow_id: Optional[int]
    overwrite_input: bool = False
