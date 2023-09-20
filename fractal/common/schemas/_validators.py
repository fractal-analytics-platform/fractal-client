import os


def valstr(attribute: str, accept_none: bool = False):
    """
    Check that a string attribute is not an empty string, and remove the
    leading and trailing whitespace characters.

    If `accept_none`, the validator also accepts `None`.
    """

    def val(string: str):
        if string is None:
            if accept_none:
                return string
            else:
                raise ValueError(
                    f"String attribute '{attribute}' cannot be None"
                )
        s = string.strip()
        if not s:
            raise ValueError(f"String attribute '{attribute}' cannot be empty")
        return s

    return val


def valint(attribute: str, min_val: int = 1):
    """
    Check that an integer attribute (e.g. if it is meant to be the ID of a
    database entry) is greater or equal to min_val.
    """

    def val(integer: int):
        if integer is None:
            raise ValueError(f"Integer attribute '{attribute}' cannot be None")
        if integer < min_val:
            raise ValueError(
                f"Integer attribute '{attribute}' cannot be less than "
                f"{min_val} (given {integer})"
            )
        return integer

    return val


def val_absolute_path(attribute: str):
    """
    Check that a string attribute is an absolute path
    """

    def val(string: str):
        if string is None:
            raise ValueError(f"String attribute '{attribute}' cannot be None")
        s = string.strip()
        if not s:
            raise ValueError(f"String attribute '{attribute}' cannot be empty")
        if not os.path.isabs(s):
            raise ValueError(
                f"String attribute '{attribute}' must be an absolute path "
                f"(given '{s}')."
            )
        return s

    return val
