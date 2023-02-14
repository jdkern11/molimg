from pathlib import Path
import pytest

from molimg import helpers


def test_remove_extension():
    data = {
        "data.py": "data",
        "data/img.png": "data/img",
        str(Path("folder") / "test" / "data.png"): str(
            Path("folder") / "test" / "data"
        ),
    }
    for filename, expected_return in data.items():
        assert helpers.remove_extensions(filename) == expected_return
