import pytest
from pathlib import Path

_TEST_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = _TEST_ROOT / "data" / "dyna"


@pytest.fixture
def pdb_path():
    p = DATA_DIR / "WT_apo_CA_ChainA.pdb"
    if not p.is_file():
        pytest.skip(f"Missing test data: {p}")
    return str(p)
