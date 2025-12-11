from md_intrinsic_dimension import intrinsic_dimension
from moleculekit.molecule import Molecule
from moleculekit.projections.metriccoordinate import MetricCoordinate
import numpy as np
import pytest
from tests.conftest import TOPO_PATH, TRAJ_PATH, REF_PATH

ATOL = 0.1


@pytest.fixture
def load_mol():
    mole = Molecule(TOPO_PATH)
    mole.read(TRAJ_PATH)
    return mole


@pytest.fixture
def load_dih_local_ID():
    return np.load(REF_PATH / "local.npy")


@pytest.fixture
def load_dih_local_all_ID():
    return np.load(REF_PATH / "mean_all.npy")


@pytest.fixture
def load_dih_local_last_ID():
    return np.load(REF_PATH / "mean_last.npy")


@pytest.fixture
def load_dih_global_all_ID():
    return np.load(REF_PATH / "global_all.npy")


@pytest.fixture
def load_dih_global_last_ID():
    return np.load(REF_PATH / "global_last.npy")


@pytest.fixture
def load_coor_local_ID():
    return np.load(REF_PATH / "local_coordinates.npy")


class TestFunctionIDResults:
    def test_coordinate_metric(self, load_mol, load_coor_local_ID):
        _, _, local_id = intrinsic_dimension(
            mol=load_mol,
            projection_method="Coordinate",
            projection_kwargs={
                "atomsel": "protein and name CA",
                "refmol": Molecule(TOPO_PATH),
            },
        )
        assert np.allclose(load_coor_local_ID, local_id, atol=ATOL)

    def test_dihedrals_metric(
        self, load_mol, load_dih_local_all_ID, load_dih_local_last_ID, load_dih_local_ID
    ):
        lid, lid100, local_id = intrinsic_dimension(
            mol=load_mol, projection_method="Dihedrals"
        )
        assert np.allclose(load_dih_local_ID, local_id, atol=ATOL)
        assert np.allclose(load_dih_local_all_ID, lid, atol=ATOL)
        assert np.allclose(load_dih_local_last_ID, lid100, atol=ATOL)

    def test_global_id(self, load_mol, load_dih_global_all_ID, load_dih_global_last_ID):
        gid, gid100 = intrinsic_dimension(
            mol=load_mol, projection_method="Dihedrals", id_method="global"
        )
        assert np.allclose(load_dih_global_all_ID, gid, atol=ATOL)
        assert np.allclose(load_dih_global_last_ID, gid100, atol=ATOL)


class TestProteinLoad:
    def test_load_mol(self, load_dih_local_ID):
        mol = Molecule(TOPO_PATH)
        mol.read(TRAJ_PATH)
        _, _, local_id = intrinsic_dimension(
            mol=mol, projection_method="Dihedrals"
        )  # default local and #distances
        assert np.allclose(load_dih_local_ID, local_id, atol=ATOL)

    def test_load_topo_traj(self, load_dih_local_ID):
        _, _, local_id = intrinsic_dimension(
            topology=TOPO_PATH, trajectory=TRAJ_PATH, projection_method="Dihedrals"
        )  # default local and #distances
        assert np.allclose(load_dih_local_ID, local_id, atol=ATOL)

    def test_load_mol_topo_traj(self, load_mol, load_dih_local_ID):
        mol = Molecule(TOPO_PATH)
        mol.read(TRAJ_PATH)
        _, _, local_id = intrinsic_dimension(
            topology=TOPO_PATH,
            trajectory=TRAJ_PATH,
            mol=load_mol,
            projection_method="Dihedrals",
        )
        assert np.allclose(load_dih_local_ID, local_id, atol=ATOL)

    def test_load_missing_traj(
        self,
    ):
        with pytest.raises(FileNotFoundError, match="Trajectory file not found"):
            intrinsic_dimension(topology=TOPO_PATH, projection_method="Dihedrals")

    def test_load_missing_topo(
        self,
    ):
        with pytest.raises(FileNotFoundError, match="Topology file not found"):
            intrinsic_dimension(trajectory=TRAJ_PATH, projection_method="Dihedrals")

    def test_missing_arguments(self):
        with pytest.raises(
            FileNotFoundError, match="file not found: None"
        ):  # stops topology file following code order
            intrinsic_dimension(projection_method="Dihedrals")

    def test_wrong_path(self):
        with pytest.raises(FileNotFoundError, match="was not found."):
            intrinsic_dimension(
                topology="wrong_path_to/topology_file",
                trajectory=TRAJ_PATH,
                projection_method="Dihedrals",
            )


class TestProjections:
    
    def test_lowercase(self, load_mol):
        with pytest.raises(
            ImportError, match="Failed to import or use custom projection class"
        ):
            intrinsic_dimension(mol=load_mol, projection_method="dihedrals")

    def test_wrong_metric(self, load_mol):
        with pytest.raises(
            ImportError, match="Failed to import or use custom projection class"
        ):
            intrinsic_dimension(mol=load_mol, projection_method="DoesNotExist")

    def test_Projection_input(self, load_mol, load_coor_local_ID):
        met = MetricCoordinate(
            atomsel="protein and name CA", refmol=Molecule(TOPO_PATH)
        )
        _, _, local_id = intrinsic_dimension(mol=load_mol, projection_method=met)
        assert np.allclose(load_coor_local_ID, local_id, atol=ATOL)

    def test_1d_array(self, load_mol):
        p_1d = np.arange(100)
        with pytest.raises(
            ValueError, match="Expected 2D array, got 1D array instead:"
        ):
            intrinsic_dimension(mol=load_mol, projection_method=p_1d)

    def test_50_frames_array(
        self, load_mol
    ):  # easy fixable by altering the neighbourhood of the ID method
        p_50 = np.arange(150).reshape(50, 3)
        with pytest.raises(ValueError, match="a minimum of 101 is required."):
            intrinsic_dimension(mol=load_mol, projection_method=p_50)


class TestID:
    def test_wrong_method(self, load_mol):
        with pytest.raises(TypeError, match='id_method must be "local" or "global"'):
            intrinsic_dimension(
                mol=load_mol,
                projection_method="Dihedrals",
                id_method="NotCorrectMethod",
            )
