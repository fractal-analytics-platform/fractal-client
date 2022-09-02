def collect_tasks():
    from fractal_tasks_core import __FRACTAL_MANIFEST__

    return (
        task_manifest
        for task_manifest in __FRACTAL_MANIFEST__
        if "task" in task_manifest["resource_type"]
    )
