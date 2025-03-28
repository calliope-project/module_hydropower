"""Set of standard clio tests.

DO NOT MANUALLY MODIFY THIS FILE!
It should be updated through our templating functions.
"""

import subprocess
from pathlib import Path

import pytest
from clio_tools.data_module import ModuleInterface


@pytest.fixture(scope="module")
def module_path():
    """Parent directory of the project."""
    return Path(__file__).parent.parent


def test_interface_file(module_path):
    """The interfacing file should be correct."""
    assert ModuleInterface.from_yaml(module_path / "INTERFACE.yaml")


def test_snakemake_all_failure(module_path):
    """The snakemake 'all' rule should return an error by default."""
    process = subprocess.run(
        "snakemake", shell=True, cwd=module_path, capture_output=True
    )
    assert "This workflow must be called as a snakemake module" in str(process.stderr)


def test_snakemake_integration_testing(module_path):
    """Run a light-weight test simulating someone using this module."""
    assert subprocess.run(
        "snakemake --use-conda --dry-run results/by_shape_id/hydro_dam_cf.parquet",
        shell=True,
        check=True,
        cwd=module_path / "tests/integration",
    )
