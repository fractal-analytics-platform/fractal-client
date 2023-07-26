from pathlib import Path
from typing import Iterable
from typing import Mapping

import mkdocs_gen_files  # type: ignore[import]
from mkdocs_gen_files import Nav

from fractal.parser import parser_main
import os
from typing import Any, Optional

from devtools import debug

import sys
from pathlib import Path
sys.path.append(Path(__file__).parent.as_posix())
from parser import parse_parser



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



nav = Nav()

# Parse main parser
main = parse_parser(parser_main)

# Create output_dir
output_dir = "reference"
if not os.path.isdir(output_dir):
    os.mkdir(output_dir)
if not os.path.isdir(f"{output_dir}/fractal"):
    os.mkdir(f"{output_dir}/fractal")

# Parser level 0
nav[["fractal"]] = f"index.md"
with open(f"{output_dir}/fractal/index.md", "w") as f:
    f.write(to_markdown(main, level=0))

# Parser levels 1 and 2 (commands and subcommands)
for child in main["children"]:
    # Level 1
    name = child["name"]
    folder = f"{output_dir}/fractal/{name}"
    if not os.path.isdir(folder):
        os.mkdir(folder)
    path = f"{folder}/index.md"
    nav[["fractal", name]] = f"{name}/index.md"
    with open(path, "w") as f:
        f.write(to_markdown(child, level=0))
        if "children" not in child.keys():
            continue
        # Level 2
        for grandchild in child["children"]:
            f.write(to_markdown(grandchild, level=1, parent_cmd=name))

with mkdocs_gen_files.open("reference/fractal/SUMMARY.md", "w") as nav_file:
    debug(nav.build_literate_nav())
    nav_file.writelines(nav.build_literate_nav())


