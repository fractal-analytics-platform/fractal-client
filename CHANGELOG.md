**Note**: Numbers like (\#123) point to closed Pull Requests on the fractal repository.

# 2.18.1 (unreleased)

* Dependencies:
    * Drop `fractal-server` and `psutil` dev dependencies (\#867).
* Testing:
    * Run Fractal Server and Postgres using Docker (\#867).

# 2.18.0

This version is aligned with [fractal-server 2.18.10](https://github.com/fractal-analytics-platform/fractal-server/blob/main/CHANGELOG.md#2180).

* CLI:
    * Deprecate `group update` command (\#854).
    * Deprecate `--viewer-paths` argument of `group new` (\#854).
    * Introduce `--add-project-dir` and `--remove-project-dir` arguments of `user edit` (\#854).
    * Deprecate `--new-project-dir` argument of `user edit` (\#854).
    * Deprecate `--zarr-dir` argument of `project add-dataset` (\#xxx)
    * Add `--project-dir` and `--zarr-subfolder` arguments of `project add-dataset` (\#xxx)
    * Make `dataset edit --new-name` argument required (\#849).
* Dependencies:
    * Support Python 3.14 (\#851).
* Docs:
    * Bump docs dependencies and update mkdocs config (\#850).
* Testing:
    * Fix tests for `fractal-server=2.17.1` (\#848).
    * Add `test_job_logs_wrong_content_type` (\#849).
    * Bump `coverage` version (\#849).

# 2.17.0

This version is aligned with [fractal-server 2.17.10](https://github.com/fractal-analytics-platform/fractal-server/blob/main/CHANGELOG.md#2170).

* CLI:
    * Add `resource new` and `profile new` commands (\#844).
    * Expose `profile_id` and `project_dir` as user's properties (\#844).
    * Deprecate arguments related to obsolete `user_settings` database table (\#844).
    * Deprecate arguments related to obsolete `user_oauth.username` database column (\#844).
* Testing:
    * Include `psutil` as a dev dependency, and use it to stop the `fractal-server` testing process (\#844).

# 2.10.2

* Fix wheel file (\#841).

# 2.10.1 [yanked]

> WARNING: This PyPI release was yanked due to invalid wheel file. Use 2.10.2 instead.

This version concerns development tools (\#831 and direct commits on `main`).

* Move `dev` and `docs` from poetry groups to standard dependency groups.
* Move from `poetry-core` to `setuptools` as a build backend.
* Move from `poetry` to `uv`.
* Move main package to `src/fractal_client`.
* Move from `mhausenblas/mkdocs-deploy-gh-pages` external action to `mkdocs gh-deploy`.
Drop


# 2.10.0

This version is aligned with [fractal-server 2.15.10](https://github.com/fractal-analytics-platform/fractal-server/blob/main/CHANGELOG.md#2150).

* CLI:
    * Introduce pre/post pinned packages (\#826).
* Dependencies:
    * Add support for Python 3.13 (\#814).
    * Bump `docs` and `dev` dependencies (\#816).
* Testing:
    * Use ubuntu-24 for GitHub CI (\#816).
    * Bump `poetry` and default Python versions in GitHub actions (\#816).


# 2.9.1

* Bump versions of `python-dotenv` and `packaging` (\#803).

# 2.9.0

* Add `--fractal-server` CLI option (\#801).
* Add `--token-path` CLI option and `FRACTAL_TOKEN_PATH` env variable (\#801).
* Stop caching token, in favor of making a call to the login API for each CLI command if needed (\#801).
* Drop `AuthToken` class (\#801).
* Expose `--task-type` in task-creation command (\#801).
* Remove default `FRACTAL_SERVER="http://localhost:8000"` configuration (\#801).

# 2.8.1

This version deprecates Python 3.10.

# 2.8.0

This version is aligned with [fractal-server 2.14.0](https://github.com/fractal-analytics-platform/fractal-server/blob/main/CHANGELOG.md#2140).

* Dependencies:
    * Bump `httpx` to `0.28.*` (\#783).

# 2.7.1

* Testing:
    * Adapt to fractal-server switch to Pydantic V2 (\#781).

# 2.7.0

This version is aligned with [fractal-server 2.11.0](https://github.com/fractal-analytics-platform/fractal-server/blob/main/CHANGELOG.md#2110).

* Commands:
    * Align with fractal-server 2.11 changes of data structure and API for filters (\#776).
* Package:
    * Remove `exclude` field, and move `packages` field to `tool.poetry` section, in pyproject.toml (\#771).
* Testing:
    * Improve configuration for coverage GitHub Action step (\#772).
    * Add `persist-credentials: false` to all `actions/checkout@v4` GitHub Action steps (\#773).

# 2.6.1

* Package:
    * Move to `poetry` v2 (\#770).
    * Require `Python<3.13` (\#770).
* Testing:
    * Use default Postgres service in GitHub CI (\#761).

# 2.6.0

* Align with new task-collection endpoint in `fractal-server` 2.10.0 (\#760).
* Update versions of pre-commit hooks and add precommit GitHub Action (\#757).

# 2.5.1

* Deprecate user `cache_dir` , to align with [fractal-server 2.9.2](https://github.com/fractal-analytics-platform/fractal-server/blob/main/CHANGELOG.md#292) (\#758).

# 2.5.0

* Update task-collection commands, to align with [fractal-server 2.9.0](https://github.com/fractal-analytics-platform/fractal-server/blob/main/CHANGELOG.md#290) (\#738).
* Remove (internal) obsolete `do_not_separate_logs` argument (\#738).
* Add `group {add|remove}-user` commands, and deprecate `--new-user-ids` argument from `group update` (\#748).
* Update `user whoami --viewer-paths` to call the new dedicated [server endpoint](https://github.com/fractal-analytics-platform/fractal-server/pull/2096) (\#748).
* Add `user set-groups` commands (\#753).
* Testing:
    * Align with fractal-server 2.9.0 removal of `DB_ENGINE` variable (\#743).

# 2.4.0

> WARNING: This release has a breaking change in the `project add-dataset` command.

* Move to from positional `zarr_dir` argument to optional `--zarr-dir` argument, for `project add-dataset` (\#736).
* Add support for user-settings `project_dir`, introduced in fractal-server 2.8.0 (\#736).
* Internal:
    * Update effect of `include_logs` for task-collection check command (\#730).

# 2.3.0

> WARNING: Starting from this release, Python3.9 is not supported any more.

* Align with [`fractal-server` 2.7.0](https://fractal-analytics-platform.github.io/fractal-server/changelog/#270) (\#712).
* Align with `fractal-server 2.7.1 and drop use of pip extras (\#727).
* Remove `--new-name` and `--new-version` options from `task edit` command (\#712).
* Rename `source` into `label`, for `task collect-custom` command (\#712).
* Do not refer to obsolete task attributes `source` or `owner` (\#712, \#717).
* Add `--new-ssh-settings-json` option to `fractal user edit` (\#715).
* Add `--ssh-settings-json` option to `fractal user register` (\#722).
* Add `--private` option to task-creating commands (\#717).
* Drop `task delete` command (\#717).
* Handle missing server in `fractal version` (\#724).
* Testing:
    * Run all tests against a single `fractal-server` instance (\#717).
    * Run tests in random module order, based on `pytest-randomly` (\#717).
    * Include Python3.12 in GitHub CI (\#717).

# 2.2.1

* Support new `viewer-paths` commands in `fractal group` commands (\#709).
* Make `--new-user-ids` optional `fractal group` command (\#709).
* Add `--viewer-paths` argument in `fractal user whoami` (\#709).

# 2.2.0

* Align with [`fractal-server` 2.6.0](https://fractal-analytics-platform.github.io/fractal-server/changelog/#260) (\#705).

# 2.1.0

* Align with [`fractal-server` 2.4.0](https://fractal-analytics-platform.github.io/fractal-server/changelog/#240) (\#695).
* Add `fractal group` command (\#695).
* Testing
    * Update GitHub actions for upload/download/coverage (\#690, \#691).
    * Switch from SQLite to Postgres in CI (\#702).

# 2.0.3

* Improve `workflow import` command (\#686).

# 2.0.2

* Improve error handling in `task collect-custom` command (\#680).
* Documentation
    * Bump `mkdocstrings[python]` to 0.25.2 (\#682).

# 2.0.1

* Add new command `task collect-custom` (\#667).
* Update `poetry` version for development to 1.8.2.
* Testing:
    * Update CI for fractal-server 2.1.0 compatibility (\#655).
    * Remove obsolete folders from `tests/data` (\#656).

# 2.0.0

Major version to align with `fractal-server` API v2.

# 1.4.4

* Require user's verification to be specified when editing user's email (\#620).

# 1.4.3

* Make `fractal-client` a fully synchronous client, by removing all `async`/`await` (\#592).
* Improve handling of `AuthenticationError` and `ConnectionError`, and handle uncaught `Exception` (\#587).
* Deprecate environment variable `FRACTAL_LOGGING_LEVEL`, remove flag `--verbose` and replace it with `--debug`, improve debugging of http requests (\#597).
* Testing:
    * Adapt `job_factory` and tests to the presence of new `timestamp_created` attributes in `fractal-server` (\#589).
    * Align with `fractal-server` 1.4.3a2 (\#598).
* Documentation:
    * Add info about server/client version compatibility (\#591).
* Dependencies:
    * Update python-dotenv to `^0.21.0` (\#589).
    * Introduce automatic updates of `poetry.lock` (\#609 and commits to `main`).

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
