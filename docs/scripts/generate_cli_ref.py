import subprocess  # nosec
from pathlib import Path

CMD0 = "fractal"


def log_help(cmd: list[str]) -> str:
    out = subprocess.check_output(
        cmd + ["--help"],
        encoding="utf-8",
    )  # nosec
    return out.strip()


def get_subcommends_from_help_message(help_msg: str) -> list[str]:
    try:
        return help_msg.split("{")[1].split("}")[0].split(",")
    except IndexError:
        return []


output_dir = Path(__file__).parents[1] / "cli_reference"
output_dir.mkdir(exist_ok=True, parents=True)

with (output_dir / "index.md").open("w") as f:
    f.write("# CLI Reference\n\n")
    f.write(
        "This page shows the help screens for the `fractal` command "
        "and its subcommands.\n\n"
    )

    f.write(f"## `{CMD0}`\n\n")
    main_help = log_help([CMD0])
    f.write("```\n")
    f.write(main_help)
    f.write("\n```\n\n")

for cmd1 in get_subcommends_from_help_message(main_help):
    with (output_dir / f"{cmd1}.md").open("w") as f:
        print(cmd1)

        f.write(f"# `{CMD0} {cmd1}`\n\n")
        help_text = log_help([CMD0, cmd1])
        f.write("```\n")
        f.write(help_text)
        f.write("\n```\n\n")

        for cmd2 in get_subcommends_from_help_message(help_text):
            f.write(f"## `{CMD0} {cmd1} {cmd2}`\n\n")
            help_text = log_help([CMD0, cmd1, cmd2])
            f.write("```\n")
            f.write(help_text)
            f.write("\n```\n\n")
