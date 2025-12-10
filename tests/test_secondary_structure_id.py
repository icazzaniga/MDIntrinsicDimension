from md_intrinsic_dimension import secondary_structure_id
from moleculekit.molecule import Molecule
import numpy as np
import pandas as pd
import pytest
from tests.conftest import TOPO_PATH, TRAJ_PATH, REF_PATH


@pytest.fixture
def load_mol():  # once established, avoid multiple loadings
    mole = Molecule(TOPO_PATH)
    mole.read(TRAJ_PATH)
    return mole


@pytest.fixture
def load_mol_ref():  # once established, avoid multiple loadings
    mole = Molecule(TOPO_PATH)
    return mole


@pytest.fixture
def load_secondary_structure_ID():
    return pd.read_pickle(REF_PATH / "secondary_structure_id.pkl")


@pytest.fixture
def load_secondary_structure_ID_table():
    return pd.read_pickle(REF_PATH / "secondary_structure_id_table.pkl")


def test_control(
    load_mol,
    load_mol_ref,
    load_secondary_structure_ID,
    load_secondary_structure_ID_table,
):
    structures, tables = secondary_structure_id(
        mol=load_mol,
        mol_ref=load_mol_ref,
        projection_method="Dihedrals",
        id_method="global",
    )
    pd.testing.assert_frame_equal(
        load_secondary_structure_ID, structures, rtol=1e-5, atol=1e-8
    )
    pd.testing.assert_frame_equal(load_secondary_structure_ID_table, tables)


class TestProteinLoad:
    def test_load_mol(self, load_mol, load_mol_ref, load_secondary_structure_ID):
        mol = Molecule(TOPO_PATH)
        mol.read(TRAJ_PATH)
        structures, _ = secondary_structure_id(
            mol=load_mol,
            mol_ref=load_mol_ref,
            projection_method="Dihedrals",
            id_method="global",
        )
        pd.testing.assert_frame_equal(
            load_secondary_structure_ID, structures, rtol=1e-5, atol=1e-8
        )

    def test_load_topo_traj(self, load_mol_ref, load_secondary_structure_ID):
        structures, _ = secondary_structure_id(
            topology=TOPO_PATH,
            trajectory=TRAJ_PATH,
            mol_ref=load_mol_ref,
            projection_method="Dihedrals",
            id_method="global",
        )
        pd.testing.assert_frame_equal(
            load_secondary_structure_ID, structures, rtol=1e-5, atol=1e-8
        )

    def test_load_mol_topo_traj(
        self, load_mol, load_mol_ref, load_secondary_structure_ID
    ):  # da rivedere
        mol = Molecule(TOPO_PATH)
        mol.read(TRAJ_PATH)
        structures, _ = secondary_structure_id(
            topology=TOPO_PATH,
            trajectory=TRAJ_PATH,
            mol=load_mol,
            mol_ref=load_mol_ref,
            projection_method="Dihedrals",
            id_method="global",
        )
        pd.testing.assert_frame_equal(
            load_secondary_structure_ID, structures, rtol=1e-5, atol=1e-8
        )

    def test_load_missing_traj(self, load_mol_ref):
        with pytest.raises(FileNotFoundError, match="Trajectory file not found"):
            secondary_structure_id(
                topology=TOPO_PATH, mol_ref=load_mol_ref, projection_method="Dihedrals"
            )

    def test_load_missing_topo(self, load_mol_ref):
        with pytest.raises(FileNotFoundError, match="Topology file not found"):
            secondary_structure_id(
                trajectory=TRAJ_PATH,
                mol_ref=load_mol_ref,
                projection_method="Dihedrals",
            )

    def test_missing_arguments(self):
        with pytest.raises(
            FileNotFoundError, match="file not found: None"
        ):  # stops topology file following code order
            secondary_structure_id(projection_method="Dihedrals")


class TestSecondaryStructure:
    def test_load_missing_refmol(self, load_mol):
        with pytest.raises(
            FileNotFoundError, match="Missing reference structure for DSSP computation"
        ):
            secondary_structure_id(mol=load_mol, projection_method="Dihedrals")

    def test_load_long_refmol(self, load_mol):
        with pytest.raises(ValueError, match="ref_mol must be 1 frame long."):
            secondary_structure_id(
                mol=load_mol, mol_ref=load_mol, projection_method="Dihedrals"
            )
