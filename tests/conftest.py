import pathlib as plib
import pytest
from hplc_data_analysis.hplc import Project


# test minimal_case
test_dir: plib.Path = plib.Path(__file__).parent
minimal_case_dir = test_dir / "data_minimal_case"


@pytest.fixture
def hplcproject() -> Project:
    proj = Project(folder_path=minimal_case_dir)
    return proj
