import argparse as ap

import pytest

from fractal_client.auth._cmd_handler import get_cmd_handler


def test_get_cmd_handler(monkeypatch):
    with pytest.raises(
        SystemExit,
        match="Help message",
    ):
        get_cmd_handler(
            args=ap.Namespace(cmd=None),
            parser_help="Help message",
        )

    with pytest.raises(
        SystemExit,
        match="internal error",
    ):
        get_cmd_handler(
            args=ap.Namespace(cmd="missing-command"),
            parser_help="Help message",
        )

    cmd_handler = get_cmd_handler(
        args=ap.Namespace(cmd="version"),
        parser_help="Help message",
    )
    assert cmd_handler is not None
