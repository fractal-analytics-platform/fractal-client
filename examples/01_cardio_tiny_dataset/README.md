## First time setup

This example uses both the client (from this repo, callable from a terminal as `fractal`) and the server (from the `fractal-server` repo, callable from the terminal as `server` - but possibly renamed later). They should both be installed in the current environment. If they are not (i.e. if one of these two commands returns a `command not found` error in the terminal), you should install the missing package(s).

To download the example dataset from Zenodo, you can use (from this folder):
```bash
./fetch_test_data_from_zenodo.sh
```
Note that this currently requires the `zenodo-get` package. At the moment this is an optional dependency, so just use `pip install zenodo-get` if you don't have it available.
TODO: this will be fixed later, either by adding it as a mandatory (dev) dependency, or by bypassing this external package.

After the setup (installing missing package(s) and downloading data) repeat the steps below for every run of the example,.
Notice that restarting the server at each run is (obviously) not going to be necessary, in the future, but at the moment it is the easiest way to have the example run consistently.

## From the `server` folder

1. If needed, run `. kill_multiprocessing.sh`;
2. Run `./server_config_and_start.sh`.
Errors like "Address already in use" point towards the need of `kill_multiprocessing.sh`).
You should see output similar to
```
[...]
INFO:     Uvicorn running on http://127.0.0.1:8000 (Press CTRL+C to quit)
[...]
INFO:     Waiting for application startup.
[...]
INFO:     Application startup complete.
```

## From the `client` folder (in a new terminal)

**WARNING**: The current scripts always delete the output folder, before starting. Make sure you change this behavior when running long examples.

Run `. define_and_apply_workflow.sh`.
A few seconds after this scripts ends, there should be two zarr files in the `tmp-proj/output-ds` folder.


## Ex-post information

1. The output images can be visualized for instance via
```bash
napari --plugin napari-ome-zarr -vvv client/tmp-proj/output-ds/20200812-CardiomyocyteDifferentiation14-Cycle1.zarr
```
2. `view_info.sh` is useful (after the `define_and_apply_workflow.sh` script has started, so that the user is registered and the environment variables are set) for a list of client-based inspection tools.
3. You can use `parsl-visualize -d` from the `server` folder.
