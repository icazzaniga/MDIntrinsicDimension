from md_intrinsic_dimension import secondary_structure_id
from moleculekit.molecule import Molecule
import numpy as np
import pandas as pd
import pytest

@pytest.fixture
def load_mol():  #once established, avoid multiple loadings
    mole = Molecule('data/2f4k.pdb')
    mole.read('data/short_2f4k_f0.xtc')
    return mole 
#mean_all = np.load('test_outputs/villin_mean_all.npy')
#mean_last = np.load('test_outputs/villin_mean_last.npy')

@pytest.fixture
def load_mol_ref():  #once established, avoid multiple loadings
    mole = Molecule('data/2f4k.pdb')
    return mole 

@pytest.fixture
def load_secondary_structure_ID(): 
    return pd.read_csv('test_outputs/secondary_structure_id.csv').drop('instantaneous', axis = 1).drop('window', axis = 1)

@pytest.fixture
def load_secondary_structure_ID_table(): 
    return pd.read_csv('test_outputs/secondary_structure_id_table.csv')









def test_control(load_mol, load_mol_ref, load_secondary_structure_ID):
    structures,_ = secondary_structure_id(mol=load_mol,mol_ref=load_mol_ref, projection_method ='Dihedrals')
    structures = structures.drop('instantaneous', axis = 1).drop('window', axis = 1) 
    pd.testing.assert_frame_equal(load_secondary_structure_ID, structures, rtol=1e-5, atol=1e-8)

    #andrebbe controllato anche df instantaneous.
    #andrebbe controllato load_secondary_structure_ID_table


class TestProteinImport:
    '''
    Check correct importing setup and conditions, before slicing the protein
    '''
    def test_import_mol(self,load_mol, load_mol_ref, load_secondary_structure_ID):
                    mol = Molecule('data/2f4k.pdb')   
                    mol.read('data/short_2f4k_f0.xtc') 
                    structures, _ = secondary_structure_id(mol=load_mol, mol_ref = load_mol_ref,projection_method='Dihedrals')
                    structures=structures.drop('instantaneous', axis = 1).drop('window', axis = 1)
                    pd.testing.assert_frame_equal(load_secondary_structure_ID, structures, rtol=1e-5, atol=1e-8)

    def test_import_topo_traj(self,load_mol_ref, load_secondary_structure_ID):
                    structures,_ = secondary_structure_id(topology='data/2f4k.pdb', trajectory='data/short_2f4k_f0.xtc',mol_ref = load_mol_ref, projection_method='Dihedrals')
                    structures=structures.drop('instantaneous', axis = 1).drop('window', axis = 1)
                    pd.testing.assert_frame_equal(load_secondary_structure_ID, structures, rtol=1e-5, atol=1e-8)
            
    def test_import_mol_topo_traj(self,load_mol, load_mol_ref, load_secondary_structure_ID): #da rivedere
                    mol = Molecule('data/2f4k.pdb')   
                    mol.read('data/short_2f4k_f0.xtc')         
                    structures,_ = secondary_structure_id(topology='data/2f4k.pdb', trajectory='data/short_2f4k_f0.xtc', mol=load_mol, mol_ref = load_mol_ref,projection_method='Dihedrals')
                    structures=structures.drop('instantaneous', axis = 1).drop('window', axis = 1)
                    pd.testing.assert_frame_equal(load_secondary_structure_ID, structures, rtol=1e-5, atol=1e-8)


    def test_import_missing_traj(self, load_mol_ref): 
        with pytest.raises(TypeError, match="path should be string"):
            secondary_structure_id(topology='data/2f4k.pdb',mol_ref = load_mol_ref, projection_method='Dihedrals')

    def test_import_missing_topo(self, load_mol_ref): 
        with pytest.raises(TypeError, match="path should be string"):
            secondary_structure_id(trajectory='data/short_2f4k_f0.xtc', mol_ref = load_mol_ref,projection_method='Dihedrals')


''' #nel codice cerco solo se c'è
    def test_import_missing_traj(self, load_mol_ref): 
        with pytest.raises(FileNotFoundError, match='Trajectory file not found'):
            secondary_structure_id(topology='data/2f4k.pdb',mol_ref = load_mol_ref, projection_method='Dihedrals')

    def test_import_missing_topo(self, load_mol_ref): 
        with pytest.raises(FileNotFoundError, match='Topology file not found'):
            secondary_structure_id(trajectory='data/short_2f4k_f0.xtc', mol_ref = load_mol_ref,projection_method='Dihedrals')
'''

'''
    def test_import_none(self, load_mol_ref):   #non è nel codice
        with pytest.raises(ImportError, match=''):
            secondary_structure_id(mol_ref = load_mol_ref,projection_method='Dihedrals')
'''

class TestSecondaryStructure:                     
    def test_import_missing_refmol(self, load_mol):
        with pytest.raises(FileNotFoundError, match='Missing reference structure for DSSP computation'):
            secondary_structure_id(mol=load_mol,projection_method='Dihedrals')

    def test_import_long_refmol(self, load_mol):
        with pytest.raises(ValueError, match='ref_mol must be 1 frame long.'):
            secondary_structure_id(mol=load_mol, mol_ref=load_mol,projection_method='Dihedrals')
    

    
#aggiungere test se il numero di atomi è diverso mol e refmol.

           


   