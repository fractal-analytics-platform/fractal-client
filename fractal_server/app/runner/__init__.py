from ._common import auto_output_dataset  # noqa: F401
from ._common import validate_workflow_compatibility  # noqa: F401


RUNNER_BACKEND = "PARSL"


if RUNNER_BACKEND == "PARSL":
    from .parsl import submit_workflow
elif RUNNER_BACKEND == "process":
    from .process import submit_workflow
else:

    def no_function(*args, **kwarsg):
        raise NotImplementedError(
            f"Runner backend {RUNNER_BACKEND} not implemented"
        )

    submit_workflow = no_function
