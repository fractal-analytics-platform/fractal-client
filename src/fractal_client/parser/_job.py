import argparse as ap


def add_job_parser(subparsers):
    job_parser = subparsers.add_parser(
        "job",
        description="Job commands.",
        allow_abbrev=False,
    )
    job_subparsers = job_parser.add_subparsers(
        title="Valid sub-commands", dest="subcmd", required=True
    )

    # job list
    job_list_parser = job_subparsers.add_parser(
        "list",
        description="List jobs for given project.",
        allow_abbrev=False,
    )
    job_list_parser.add_argument(
        "project_id",
        type=int,
        help="Project ID.",
    )

    # job show
    job_show_parser = job_subparsers.add_parser(
        "show",
        description="Query status of workflow-execution job.",
        argument_default=ap.SUPPRESS,
        allow_abbrev=False,
    )
    job_show_parser.add_argument("project_id", type=int, help="Project ID.")
    job_show_parser.add_argument("job_id", type=int, help="Job ID.")

    # job download-logs
    job_download_logs_parser = job_subparsers.add_parser(
        "download-logs",
        description="Download full folder of workflow-execution job.",
        allow_abbrev=False,
    )
    job_download_logs_parser.add_argument(
        "project_id", type=int, help="Project ID."
    )
    job_download_logs_parser.add_argument("job_id", type=int, help="Job ID.")
    job_download_logs_parser.add_argument(
        "--output",
        dest="output_folder",
        help="Path of the output folder.",
        required=True,
    )

    # job stop
    job_stop_parser = job_subparsers.add_parser(
        "stop",
        description="Stop workflow-execution job.",
        allow_abbrev=False,
    )
    job_stop_parser.add_argument("project_id", type=int, help="Project ID.")
    job_stop_parser.add_argument("job_id", type=int, help="Job ID.")

    # job submit
    job_submit_parser = job_subparsers.add_parser(
        "submit",
        description="Submit a job.",
        argument_default=ap.SUPPRESS,
        allow_abbrev=False,
    )

    job_submit_parser.add_argument("project_id", type=int)
    job_submit_parser.add_argument("workflow_id", type=int)
    job_submit_parser.add_argument("dataset_id", type=int)
    job_submit_parser.add_argument(
        "--start",
        dest="first_task_index",
        type=int,
        help=(
            "Positional index of the first task to be executed"
            " (starting from 0)."
        ),
        required=False,
    )
    job_submit_parser.add_argument(
        "--end",
        dest="last_task_index",
        type=int,
        help=(
            "Positional index of the last task to be executed"
            " (starting from 0)."
        ),
        required=False,
    )
    job_submit_parser.add_argument(
        "-w",
        "--worker-init",
        help="Command to be run before starting a worker.",
    )
    job_submit_parser.add_argument(
        "--attribute-filters-json",
        help=(
            "Path to JSON file with the attribute filters for this job submission."
        ),
        required=False,
    )
    job_submit_parser.add_argument(
        "--type-filters-json",
        help=(
            "Path to JSON file with the type filters for this job submission."
        ),
        required=False,
    )
