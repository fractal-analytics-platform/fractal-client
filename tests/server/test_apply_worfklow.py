from devtools import debug


async def test_apply_workflow(
    client, MockCurrentUser, project_factory, dataset_factory, resource_factory
):
    """
    GIVEN a dataset and a non-trivial workflow
    WHEN one applys the workflow to the dataset
    THEN the job is correctly submitted and executed
    """
    async with MockCurrentUser(persist=True) as user:
        prj = await project_factory(user)
        ds = await dataset_factory(prj)

        resource = await resource_factory(ds)

        debug(ds)
        debug(resource)
