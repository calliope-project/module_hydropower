"""A collection of heavier tests to run locally.

Should not be run during integration tests.
Only use to solve bugs or test new features!
"""

import shutil
import subprocess
from pathlib import Path

import pytest


def setup_workflow(workflow_path, dest_path, case):
    """Copy workflow-relevant folders and files to a temporary location."""
    shutil.copytree(workflow_path / "workflow", dest_path / "workflow")
    files_to_copy = [
        (f"tests/files/{case}/config.yaml", "config/config.yaml"),
        (
            f"tests/files/{case}/powerplants.parquet",
            "resources/user/powerplants.parquet",
        ),
        (f"tests/files/{case}/shapes.parquet", "resources/user/shapes.parquet"),
        (
            "tests/files/national_generation.parquet",
            "resources/user/national_generation.parquet",
        ),
    ]
    for src, dest in files_to_copy:
        dest_file = Path(dest_path / dest)
        dest_file.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(workflow_path / src, dest_file)


@pytest.fixture(scope="module")
def module_path():
    """Parent directory of the project."""
    return Path(__file__).parent.parent


def test_assorted_europe(module_path, tmp_path):
    """Test a small subset of european nations.

    A buffer of 10km was used to select powerplants.
    Thus, a few powerplants outside the boundary should be adjusted, and a few dropped.
    """
    setup_workflow(module_path, tmp_path, "assorted_europe")
    result_file = "results/adjusted_powerplants.parquet"
    subprocess.run(
        f"snakemake --cores 4 {result_file}", shell=True, check=True, cwd=tmp_path
    )
    assert Path(tmp_path / result_file).exists()
