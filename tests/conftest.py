from pathlib import Path

TEST_DATA_DIR = Path(__file__).parent / "data"

TOPO_PATH = str(TEST_DATA_DIR / "2hbaA00.pdb")
TRAJ_PATH = str(TEST_DATA_DIR / "2hbaA00_320_0.xtc")

REF_PATH = Path(__file__).parent / "test_outputs"
