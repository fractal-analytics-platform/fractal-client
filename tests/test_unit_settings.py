import pytest
from devtools import debug
from pydantic import ValidationError

from fractal_server.config import Settings
from fractal_server.dependency_injection import Inject


def test_settings_injection(override_settings):
    """
    GIVEN an Inject object with a Settigns object registered to it
    WHEN I ask for the Settings object
    THEN it gets returned
    """
    settings = Inject(Settings)
    debug(settings)
    assert isinstance(settings, Settings)


@pytest.mark.parametrize(
    ("settings", "raises"),
    [
        (Settings(), True),
        (
            Settings(
                DEPLOYMENT_TYPE="development",
                JWT_SECRET_KEY="secret",
                SQLITE_PATH="path",
            ),
            False,
        ),
    ],
)
def test_settings_check(settings: Settings, raises: bool):
    debug(settings)
    if raises:
        with pytest.raises(ValidationError):
            settings.check()
    else:
        settings.check()
