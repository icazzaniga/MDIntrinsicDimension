from md_intrinsic_dimension import section_id
from moleculekit.molecule import Molecule
import numpy as np
import pandas as pd
import pytest
from tests.conftest import TOPO_PATH, TRAJ_PATH, REF_PATH


@pytest.fixture
def load_mol():
    mole = Molecule(TOPO_PATH)
    mole.read(TRAJ_PATH)
    return mole


@pytest.fixture
def load_section_ID():
    return pd.read_pickle(REF_PATH / "section_id.pkl")


def test_control(
    load_mol, load_section_ID
):  # if all is imput correctly, the call works
    sections = section_id(
        mol=load_mol, projection_method="Dihedrals", id_method="global"
    )
    pd.testing.assert_frame_equal(load_section_ID, sections, rtol=1e-5, atol=1e-8)


class TestProteinLoad:
    def test_load_mol(self, load_section_ID):
        mol = Molecule(TOPO_PATH)
        mol.read(TRAJ_PATH)
        sections = section_id(
            mol=mol, projection_method="Dihedrals", id_method="global"
        )
        pd.testing.assert_frame_equal(load_section_ID, sections, rtol=1e-5, atol=1e-8)

    def test_load_topo_traj(self, load_section_ID):
        sections = section_id(
            topology=TOPO_PATH,
            trajectory=TRAJ_PATH,
            projection_method="Dihedrals",
            id_method="global",
        )
        pd.testing.assert_frame_equal(load_section_ID, sections, rtol=1e-5, atol=1e-8)

    def test_load_mol_topo_traj(self, load_mol, load_section_ID):  # da rivedere
        mol = Molecule(TOPO_PATH)
        mol.read(TRAJ_PATH)
        sections = section_id(
            topology=TOPO_PATH,
            trajectory=TRAJ_PATH,
            mol=load_mol,
            projection_method="Dihedrals",
            id_method="global",
        )
        pd.testing.assert_frame_equal(load_section_ID, sections, rtol=1e-5, atol=1e-8)

    def test_load_missing_traj(self):
        with pytest.raises(FileNotFoundError, match="Trajectory file not found"):
            section_id(topology=TOPO_PATH, projection_method="Dihedrals")

    def test_load_missing_topo(self):
        with pytest.raises(FileNotFoundError, match="Topology file not found"):
            section_id(trajectory=TRAJ_PATH, projection_method="Dihedrals")

    def test_missing_arguments(self):
        with pytest.raises(FileNotFoundError, match="file not found: None"):
            section_id(projection_method="Dihedrals")


class TestSections:
    def test_short_window(self, load_mol):
        with pytest.raises(ValueError, match="`window_size` must be > 1."):
            section_id(mol=load_mol, projection_method="Dihedrals", window_size=1)

    def test_wrong_method(self, load_mol):
        with pytest.raises(TypeError, match='id_method must be "local" or "global"'):
            section_id(mol=load_mol, projection_method="Dihedrals", id_method="Local")
