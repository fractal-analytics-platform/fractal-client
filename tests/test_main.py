import pytest

from fractal_client.client import main
from fractal_client.interface import Interface


def test_unit_main(monkeypatch):
    """
    Run the `main` function.

    NOTE: Mocking `handle()` is necessary because there is no
    appropriate `sys.argv`.
    """
    import fractal_client.client

    monkeypatch.setattr(
        fractal_client.client,
        "handle",
        lambda: Interface(data="data", retcode=0),
    )
    with pytest.raises(SystemExit):
        main()
