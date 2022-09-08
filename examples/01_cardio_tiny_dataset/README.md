The first time you want to run this example, use:
```
git fetch
git checkout dev
git pull
poetry install
```
and also (from this folder):
```bash
./fetch_test_data_from_zenodo.sh
```

Then, for every run of the example, repeat the steps below.
Notice that restarting the server at each run is (obviously) not going to be necessary, in the future, but at the moment it is the easiest way to have the example run.

### From the `server` folder

1. If needed, run `. kill_multiprocessing.sh`;
2. Run `. server_setup.sh`;
3. Run `. server_start.sh` (message errors like "Address already in use" point towards the need of `kill_multiprocessing.sh`). You should see output similar to
```
[...]
INFO:     Uvicorn running on http://127.0.0.1:8000 (Press CTRL+C to quit)
[...]
INFO:     Waiting for application startup.
[...]
INFO:     Application startup complete.
```

### From the `client` folder (in a new terminal)

1. Run `. define_and_apply_workflow.sh`. A few seconds after this scripts ends, there should be two zarr files in the `tmp-proj/output-ds` folder.
2. `view_info.sh` is useful (after the first script started, so that the user is registered and the environment variables are set) for a list of client-based inspection tools.
3. The output images can be visualized for instance via
```bash
napari --plugin napari-ome-zarr -vvv client/tmp-proj/output-ds/20200812-CardiomyocyteDifferentiation14-Cycle1.zarr
```
