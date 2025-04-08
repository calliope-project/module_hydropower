"""A collection of heavier tests to run locally.

Useful for debugging or testing new features.

Important things to consider:
 - These tests should be run individually to avoid excessive workloads.
 - Should only be run locally, as they are likely too heavy for Github's CI.
"""

import shutil
import subprocess
from pathlib import Path

import pytest


def setup_workflow(workflow_path: Path, dest_path: Path, case: str):
    """Copy workflow-relevant folders and files to a temporary location."""
    if dest_path.exists():
        shutil.rmtree(dest_path)
    shutil.copytree(workflow_path / "workflow", dest_path / "workflow")
    files_to_copy = [
        ("config/config.yaml", "config/config.yaml"),
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


def test_assorted_europe(module_path):
    """Test a small subset of european nations.

    A buffer of 10km was used to select powerplants.
    Thus, a few powerplants outside the boundary should be adjusted, and a few dropped.
    """
    case = "assorted_europe"
    tmp_path = Path(f"temp/{case}")
    setup_workflow(module_path, tmp_path, case)
    result_file = "results/by_shape_id/hydro_dam_cf.parquet"
    subprocess.run(
        f"snakemake --cores 4 {result_file}", shell=True, check=True, cwd=tmp_path
    )
    assert Path(tmp_path / result_file).exists()
