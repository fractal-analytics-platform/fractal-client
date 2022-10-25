from typing import Any
from typing import Callable
from typing import Dict
from typing import TypeVar
from typing import Union

T = TypeVar("T")
TT = TypeVar("TT")


_instance_count = 0
_UNDEFINED = "__UNDEFINED__"


class InjectionError(KeyError):
    pass


class _Inject:
    _dependencies: Dict[Any, Any] = {}

    def __init__(self):
        global _instance_count
        if _instance_count == 1:
            raise RuntimeError("You must only instance this class once")
        _instance_count += 1

    @classmethod
    def __call__(
        cls, Type: Callable[..., T], default: TT = _UNDEFINED  # type: ignore
    ) -> Union[T, TT]:
        """
        Return the dependency corresponding to Type

        If the dependency was never registered and `default` was provided,
        return the defaul value, otherwise raise an InjectionError.
        """
        try:
            return cls._dependencies[Type]
        except KeyError:
            if default != _UNDEFINED:
                return default
            else:
                raise InjectionError(f"No dependency for `{Type.__name__}`")

    @classmethod
    def pop(
        cls, Type: Callable[..., T], default: TT = _UNDEFINED  # type: ignore
    ) -> Union[T, TT]:
        """
        Pop the dependency corresponding to Type
        """
        try:
            return cls._dependencies.pop(Type)
        except KeyError:
            if default != _UNDEFINED:
                return default
            else:
                raise InjectionError(f"No dependency for `{Type.__name__}`")

    @classmethod
    def register(cls, Type: Callable[..., T], value: T) -> None:
        """
        Register a dependency
        """
        cls._dependencies[Type] = value


# NOTE: This is a singleton instance
Inject = _Inject()
