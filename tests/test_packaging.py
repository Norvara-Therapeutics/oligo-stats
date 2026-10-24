"""Invariants of the distribution itself."""

import re
from pathlib import Path

import pytest

from oligo_stats import __version__

ROOT = Path(__file__).resolve().parent.parent
PYPROJECT = ROOT / "pyproject.toml"


def _declared_version() -> str:
    """Read the version from pyproject without tomllib, which needs 3.11."""
    for line in PYPROJECT.read_text(encoding="utf-8").splitlines():
        match = re.fullmatch(r'version = "([^"]+)"', line.strip())
        if match:
            return match.group(1)
    raise AssertionError("pyproject.toml declares no version")


@pytest.mark.skipif(not PYPROJECT.is_file(), reason="not an installed tree")
def test_the_packaged_and_importable_versions_agree():
    """A drift here ships a wheel whose --version lies about what it is."""
    assert _declared_version() == __version__


@pytest.mark.skipif(not PYPROJECT.is_file(), reason="not an installed tree")
def test_the_manifest_ships_the_tests_and_examples_the_readme_cites():
    manifest = (ROOT / "MANIFEST.in").read_text(encoding="utf-8")
    assert "recursive-include tests *.py" in manifest
    for suffix in ("*.fasta", "*.tsv"):
        assert suffix in manifest
    for name in ("primers.fasta", "order_sheet.tsv"):
        assert (ROOT / "examples" / name).is_file()
