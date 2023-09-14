import pytest
from devtools import debug
from packaging import version

from fractal.cmd._aux_task_caching import _loose_version_parse


# string, major, minor, micro, is_prerelease
VERSIONS = [
    ("0.10.0c0", 0, 10, 0, True),
    ("0.10.0b4", 0, 10, 0, True),
    ("0.10.0", 0, 10, 0, False),
    ("0.10.0alpha3", 0, 10, 0, True),
    ("0.10.0a2", 0, 10, 0, True),
    ("1.0.0", 1, 0, 0, False),
    ("0.10.0a0", 0, 10, 0, True),
    ("1.0.0rc4.dev7", 1, 0, 0, True),
    ("0.10.0beta5", 0, 10, 0, True),
    ("0.10.0alpha0", 0, 10, 0, True),
    ("3.2", 3, 2, 0, False),
    ("2", 2, 0, 0, False),
]

SORTED_VERSIONS = [
    "0.10.0a0",
    "0.10.0alpha0",
    "0.10.0a2",
    "0.10.0alpha3",
    "0.10.0b4",
    "0.10.0beta5",
    "0.10.0c0",
    "0.10.0",
    "1.0.0rc4.dev7",
    "1.0.0",
    "2",
    "3.2",
]


def test_version_parsing():
    for (v_string, major, minor, micro, is_prerelease) in VERSIONS:
        v = version.parse(v_string)
        debug(v_string, v.major, v.minor, v.micro, v.is_prerelease, v.pre)
        if major is not None:
            assert v.major == major
        if minor is not None:
            assert v.minor == minor
        if micro is not None:
            assert v.micro == micro
        if is_prerelease is not None:
            assert v.is_prerelease == is_prerelease

    with pytest.raises(version.InvalidVersion):
        version.parse("invalid")


def test_version_sorting():

    sorted_versions = sorted([v[0] for v in VERSIONS], key=version.parse)
    debug(sorted_versions)
    assert sorted_versions == SORTED_VERSIONS


def test_max_with_loose_version_parse():
    versions = ["invalid_1"] + [v[0] for v in VERSIONS] + ["invalid_2"]
    sorted_versions = sorted(versions, key=_loose_version_parse)
    debug(sorted_versions)
    assert sorted_versions == ["invalid_1", "invalid_2"] + SORTED_VERSIONS
