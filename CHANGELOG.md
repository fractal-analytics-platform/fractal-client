**Note**: Numbers like (\#123) point to closed Pull Requests on the fractal repository.

# 1.4.3 (unreleased)

* Improve handling of `AuthenticationError` and `ConnectionError`, and handle uncaught `Exception` (\#587).
* Remove all `async`/`await` (\#592).
* Testing:
    * Adapt `job_factory` and tests to the presence of new `timestamp_created` attributes in `fractal-server` (\#589).
* Documentation:
    * Add info about server/client version compatibility (\#591).

# 1.4.2

* Always make new users verified, within `user register` command (\#580).
* Expose verification-related features in `user edit` command (\#580).
* Update expected status code of stop-job to 202 (\#585).
* Testing:
    * Adapt `job_factory` to new strict response-validation models in `fractal-server` (\#580).
    * Adapt `job_factory` by including the `timestamp_created` attribute (\#582).

# 1.4.1

We are skipping this patch release, to remain aligned with `fractal-server` version.

# 1.4.0

* Align with [`fractal-server` 1.4.0](https://fractal-analytics-platform.github.io/fractal-server/changelog/#140) (\#573).
* Testing:
    * Use ubuntu-22 for GitHub CI (commit e1c8bd3da670c24a0ac48b1163cd1c7833746aaf).
* Development:
    * Do not include `sqlmodel` dependency explicitly (\#577).
    * Use poetry 1.7.1 (\#577).

# 1.3.7

This release is up-to-date with `fractal-server` 1.3.12.

* Remove `project new --dataset` argument (\#566).
* Add `project new --make-read-only` argument (\#566).
* Do not use deprecated fractal-server `deployment_type` variable (\#566).
* Align tests with `fractal-server` 1.3.12, which deprecates the legacy history version (\#569).
* Automate procedure for publishing package to PyPI (\#561).

# 1.3.6

* Main features:
    * Remove client-side validation of API request/response bodies (\#551).
    * Make error messages more readable for request-body validation errors (\#551).
    * Include `--batch` option for workflow apply command (commit 06c7ff0e92602f08a98097d3597a8ce39c6ae1a8).
    * Revamp `config.py`, making `Settings` a standard Python class (\#559).
* Package and repository:
    * Rename repository from `fractal` to `fractal-client`.
    * Change package name from `fractal` to `fractal-client` (\#555).
    * Remove `fractal-common` submodule (\#551).
    * Remove `pydantic` dependency (\#559).
* Tests:
    * Review tests: make them stateless, improve handling of cache, remove obsolete fixtures (\#559).

# 1.3.5

* Implement more robust sorting of versions (e.g. in the presence of pre-releases), via `packaging.version` (\#550).

# 1.3.4

* Add new commands `dataset status` and `dataset history` (\#544).
* Align with fractal-server/fractal-common new `TaskUpdate` model, that accepts `version=None` (\#540).
* Align with fractal-server/fractal-common new attributes in the Task model (\#534).
* Align with fractal-common new `ApplyWorkflowRead` model, with optional `workflow_dump` (\#547).
* Move documentation from sphinx to mkdocs (\#532).

# 1.3.3

* Support `workflow apply --start/--end` arguments for submitting a workflow subset (\#527).
* Exclude `common/tests` and other files from build (\#523).
* Remove obsolete folders from `tests/data` (\#526).

# 1.3.2

* Fix wrong build in 1.3.1 (include missing `fractal.common` submodule).

# 1.3.1

WARNING: wrong build, do not use

* Pin Pydantic to V1 (\#520).

# 1.3.0

* Align with [fractal-server 1.3.0](https://fractal-analytics-platform.github.io/fractal-server/changelog/#130), by updating all relevant endpoint path/query/body parameters (\#479).
* Add `fractal job stop` command (\#485).
* Add `fractal task delete` command (\#510).
* Add task ID/name/version disambiguation to `task edit` and `workflow add-task` (\#499).
* Specific changes to existing commands:
    * Make `project_id` a required positional argument of `fractal {workflow,dataset,job}` commands (\#479).
    * For `edit` commands, always prepend the new arguments with `new`, e.g. as in `task edit ... --new-version` (\#498).
    * Add `--type` optional argument to `fractal dataset new` command (\#479).
    * For `fractal workflow apply`:
        * Transform `project_id` from option to positional argument (\#479).
        * Make `output_dataset_id` a required positional argument (\#483).
    * Add `--username/--new-username` to `fractal user` subcommands (\#493).
    * Remove `--private` option for `fractal task collect` (\#493).
    * Add `--version` to `fractal task {new,edit}` subcommands (\#493).
    * Split `task-id-or-name` argument of `task edit` and `workflow add-task` into two arguments (\#504).
    * Add `--pinned-dependency` argument to `task collect` (\#508).
    * Add `--args-schema` and `--args-schema-version` arguments to `task new` command (\#511).
    * Add `--new-args-schema` and `--new-args-schema-version` arguments to `task edit` command (\#511).
    * Raise warning when importing/exporting workflows with custom tasks (\#513).
* Package and repository:
    * Fix a bug in tests, by starting the fractal-server FastAPI app in a more standard way (\#481).
    * Require pydantic version to be `>=1.10.8` (\#486, \#490).
    * Make `sqlmodel` a development depedency (\#493).
    * Improve handling of a `ConnectError` in the CI (\#502).
    * Remove arbitrary `kwargs` from internal functions (\#503).
    * Align with latest fractal-server version and update tests (\#517).

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
