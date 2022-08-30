The first time you need to run this example:
```
git fetch
git checkout server-newtasks
git pull
poetry install
```

Then, for every run of the example, repeat all steps below. Notice that restarting the server at each run is (obviously) not going to be necessary, later, but for this example it works this way.

From `server_dir`:
1. If needed, run `. kill_multiprocessing.sh`;
2. Run `. server_setup.sh`;
3. Run `. server_start.sh` (message errors like "Address already in use" point towards the need of `kill_multiprocessing.sh`).

From `user_dir` (in a new terminal):
1. Run `. define_and_apply_workflow.sh`;
2. `view_info.sh` is useful (after the first script started, so that the user is registered and the environment variables are set) for a list of inspection tools.
