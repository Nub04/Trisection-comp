"""Test-suite configuration: locate the C++ binaries."""
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "python"))


def pytest_collection_modifyitems(config, items):
    pass
