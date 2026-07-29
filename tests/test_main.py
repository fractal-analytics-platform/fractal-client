import pytest

from fractal_client.interface import Interface
from fractal_client.main import main


def test_unit_main(monkeypatch):
    """
    Run the `main` function.

    NOTE: Mocking `handle()` is necessary because there is no
    appropriate `sys.argv`.
    """
    import fractal_client.main

    monkeypatch.setattr(
        fractal_client.main,
        "handle",
        lambda _: Interface(data="data", retcode=0),
    )
    with pytest.raises(SystemExit):
        main()
