def collect_tasks():
    from fractal.tasks import __FRACTAL_MANIFEST__

    return (
        task_manifest
        for task_manifest in __FRACTAL_MANIFEST__
        if task_manifest["resource_type"] == "task"
    )
