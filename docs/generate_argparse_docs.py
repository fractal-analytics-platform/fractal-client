from fractal.parser import parser_main
import os
from parser import parse_parser
from typing import Any, Optional


def to_markdown(
        data: dict[str, str],
        level: int,
        parent_cmd: Optional[str] = None,
        ) -> str:
    """
    Given a `data` object with keys `name`, `description` and `usage`, produce
    a markdown string.
    """
    name = data["name"]
    description = data["description"]
    usage = data["usage"]

    # Create MarkDown string for title
    if parent_cmd:
        title_str = "#" * (level + 1) + f" {parent_cmd.title()} {name.title()}\n"
    else:
        title_str = "#" * (level + 1) + f" {name.title()}\n"

    # Create MarkDown string for description
    description_str = f"{description}\n"

    # Create MarkDown string for usage code block
    usage_str = f"```\n{usage}\n```\n"

    # Combine strings together
    md_string = "\n".join((
            title_str,
            description_str,
            usage_str,
            )) + "\n"

    return md_string


def main(output_dir: str = "reference"):

    # Create output_dir
    if os.path.isdir(output_dir):
        raise RuntimeError(f"{output_dir=} already exists. Exit.")
    os.mkdir(output_dir)

    # Initialize summary file
    f_summary = open(f"{output_dir}/SUMMARY.md", "w")

    # Browse first two levels of parser (commands and subcommands)
    main = parse_parser(parser_main)
    for child in main["children"]:
        name = child["name"]
        f_summary.write(f"[{name}]({name}.md)\n")

        with open(f"{output_dir}/{name}.md", "w") as f:
            f.write(to_markdown(child, level=0))
            if "children" not in child.keys():
                continue
            for grandchild in child["children"]:
                f.write(to_markdown(grandchild, level=1, parent_cmd=name))

    f_summary.close()

main()
