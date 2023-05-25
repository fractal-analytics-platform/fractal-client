**Note**: Numbers like (\#123) point to closed Pull Requests on the fractal repository.

# 1.3.0

* Align with [fractal-server 1.3.0](https://fractal-analytics-platform.github.io/fractal-server/changelog/#130), by updating all relevant endpoint path/query/body parameters (\#479).
* Make `project_id` a required positional argument of `fractal {workflow,dataset,job}` commands (\#479).
* For `fractal workflow apply`, transform `project_id` from option to positional argument (\#479).
* Add `fractal job stop` command (\#485).
* Add `--type` optional argument to `fractal dataset new` command (\#479).
* Make `output_dataset_id` a required positional argument of `fractal workflow apply` (\#483).
* Fix a bug in tests, by starting the fractal-server FastAPI app in a more standard way (\#481).
* Require pydantic version to be `>=1.10.8` (\#486, \#490).
* Remove `--private` option for `fractal task collect` (\#493).
* Add `--username/--new-username` to `fractal user` subcommands (\#493).
* Add `--version` to `fractal task` subcommands (\#493).

# 1.2.0

* Align with [fractal-server 1.2.0](https://fractal-analytics-platform.github.io/fractal-server/changelog/#120) (\#472).

# 1.1.0

* Align with [fractal-server 1.1.0](https://fractal-analytics-platform.github.io/fractal-server/changelog/#110) (\#446).
* Improve validation of API request payloads (\#447).
* Drop support for python 3.8 (\#438).
* Update `_TaskBase` schema from `fractal-common` (\#431).
* Update `DatasetUpdate` schema (\#461).
* Update `fractal task edit` command (\#439 and \#461).
* Add `fractal project edit` command (\#465).
* Improve task-collection log formatting (\#443).
* Disable `argparse` abbreviation for CLI commands (\#441).

# 1.0.5

* Minor updates to `fractal workflow export` (\#429).

# 1.0.4

* Add `fractal workflow {import,export}` commands (\#426).
* Remove `--project-id` argument from `fractal workflow edit` commands (\#426).

# 1.0.3

* Add `fractal task new` command (\#421).
* Remove obsolete `-j/--json` argument from `fractal` command (\#421).

# 1.0.2

* Remove obsolete references to SLURM user, either as a CLI argument or an environment variable (\#419).

# 1.0.1

* Make `FRACTAL_USER/FRACTAL_PASSWORD` env variables optional (\#416).

# 1.0.0

* First version in this CHANGELOG.
