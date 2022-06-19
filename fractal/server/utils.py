from warnings import warn as _warn

from .config import DeploymentType
from .config import settings


def warn(message):
    """
    Make sure that warnings do not make their way to staing and production
    """
    if settings.DEPLOYMENT_TYPE in [
        DeploymentType.testing,
        DeploymentType.development,
    ]:
        _warn(message, RuntimeWarning)
    else:
        raise RuntimeError(message)
