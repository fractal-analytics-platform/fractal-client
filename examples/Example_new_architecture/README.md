The first time you want to run this example, use:
```
git fetch
git checkout server-newtasks
git pull
poetry install
```

Then, for every run of the example, repeat the steps below.
Notice that restarting the server at each run is (obviously) not going to be necessary, in the future, but at the moment it is the easiest way to have the example run.

From `server_dir`:
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

From `user_dir` (in a new terminal):
1. Run `. define_and_apply_workflow.sh`. The last item in the output should look like
```
    res.json(): {
        'status': 'submitted',
    } (dict) len=1
```
2. `view_info.sh` is useful (after the first script started, so that the user is registered and the environment variables are set) for a list of client-based inspection tools.
3. The output images can be visualized via
```bash
napari --plugin napari-ome-zarr -vvv tmp-proj/output-ds/myplate.zarr
napari --plugin napari-ome-zarr -vvv tmp-proj/output-ds/myplate_mip.zarr
```
