import sys
from pathlib import Path
from textwrap import fill
from typing import Any
from typing import Optional

import mkdocs_gen_files  # type: ignore[import]
from mkdocs_gen_files import Nav
from parser import parse_parser

from fractal.parser import parser_main

sys.path.append(Path(__file__).parent.as_posix())


# FIXME: mutually exclusive group??


def to_markdown(
    data: dict[str, Any],
    level: int,
    parent_cmd: Optional[str] = None,
) -> str:
    """
    Given a `data` object with keys `name`, `description` and `usage`, produce
    a markdown string.
    """

    # Create MarkDown string for title
    name = data["name"]
    if parent_cmd:
        title_str = "#" * (level + 2) + f" {parent_cmd} {name}\n"
    else:
        title_str = "#" * (level + 1) + f" {name}\n"

    # Create MarkDown string for description
    description = data["description"]
    description_str = f"{description}\n"

    # Create MarkDown string for usage code block
    usage = data["usage"].replace(
        "gen_ref_pages.py", "fractal"
    )  # FIXME: why is this needed??  # noqa
    while "  " in usage:
        usage = usage.replace("  ", " ")
    usage = fill(
        usage,
        width=80,
        initial_indent="",
        subsequent_indent=(" " * 16),
        break_on_hyphens=False,
    )
    usage_str = f"```\n{usage}\n```\n"

    # Create MarkDown string for action groups
    action_groups_strings = []
    if "action_groups" in data.keys():
        for group in data["action_groups"]:
            title = group["title"]

            # FIXME: how to handle these specific cases? Is it actually useful?
            if title == "Commands":
                continue
            elif title == "Valid sub-commands":
                action_groups_strings.append(
                    "#" * (level + 2) + " Sub-commands"
                )
            elif title in [
                "Named Arguments",
                "Positional Arguments",
            ]:
                options = group["options"]
                action_groups_strings.append("#" * (level + 3) + f" {title}\n")
                # FIXME: add here (?): understand whether there are
                # mutually-exclusive groups, and mark them as such
                for opt in options:
                    opt_name = ",".join(opt["name"])
                    opt_help = opt["help"]
                    # FIXME: if opt_name in some mutually-exclusive group:
                    #     opt_help += "(incompatible with SOMEOTHEROPTION)
                    action_groups_strings.append(
                        f"- **`{opt_name}`**: {opt_help}\n"
                    )
                    # FIXME: add here "required" mark?
                    # FIXME: add here "default" mark?
            else:
                raise NotImplementedError(title)

    action_groups_str = "\n".join(action_groups_strings)

    # Combine strings together
    md_string = (
        "\n".join(
            (
                title_str,
                description_str,
                usage_str,
                action_groups_str,
            )
        )
        + "\n"
    )

    return md_string


nav = Nav()

# Parse main parser
main = parse_parser(parser_main)

# Parser level 0
nav[["fractal"]] = "fractal/index.md"
main["name"] = "fractal"
with mkdocs_gen_files.open("reference/fractal/index.md", "w") as f:
    f.write(to_markdown(main, level=0))

# Parser levels 1 and 2 (commands and subcommands)
for child in main["children"]:
    # Level 1
    name = child["name"]

    nav[["fractal", name]] = f"fractal/{name}/index.md"
    with mkdocs_gen_files.open(f"reference/fractal/{name}/index.md", "w") as f:
        f.write(to_markdown(child, level=0))
        if "children" not in child.keys():
            continue
        # Level 2
        for grandchild in child["children"]:
            f.write(to_markdown(grandchild, level=1, parent_cmd=name))

with mkdocs_gen_files.open("reference/SUMMARY.md", "w") as nav_file:
    nav_file.writelines(nav.build_literate_nav())
