from sqlmodel import SQLModel


class WorkflowBase(SQLModel):
    name: str
    project_id: int


class WorkflowRead(WorkflowBase):
    id: int
