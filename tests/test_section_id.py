from md_intrinsic_dimension import section_id
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
def load_section_ID(): 
    return pd.read_csv('test_outputs/section_id.csv').drop('instantaneous', axis = 1)






def test_control(load_mol, load_section_ID): #if all is imput correctly, the call works
    sections = section_id(mol=load_mol, projection_method ='Dihedrals').drop('instantaneous', axis = 1)
    pd.testing.assert_frame_equal(load_section_ID, sections, rtol=1e-5, atol=1e-8)

    #andrebbe controllato anche df instantaneous.


class TestProteinImport:

    '''
    Check correct importing setup and conditions, before slicing the protein
    '''
  
    def test_import_mol(self, load_section_ID):
                    mol = Molecule('data/2f4k.pdb')   
                    mol.read('data/short_2f4k_f0.xtc') 
                    sections = section_id(mol=mol, projection_method='Dihedrals').drop('instantaneous', axis = 1)
                    pd.testing.assert_frame_equal(load_section_ID, sections, rtol=1e-5, atol=1e-8)

    def test_import_topo_traj(self, load_section_ID):
                    sections = section_id(topology='data/2f4k.pdb', trajectory='data/short_2f4k_f0.xtc', projection_method='Dihedrals').drop('instantaneous', axis = 1) 
                    pd.testing.assert_frame_equal(load_section_ID.drop('instantaneous', axis = 1), sections, rtol=1e-5, atol=1e-8)
            
    def test_import_mol_topo_traj(self,load_mol, load_section_ID): #da rivedere
                    mol = Molecule('data/2f4k.pdb')   
                    mol.read('data/short_2f4k_f0.xtc')         
                    sections = section_id(topology='data/2f4k.pdb', trajectory='data/short_2f4k_f0.xtc', mol=load_mol, projection_method='Dihedrals').drop('instantaneous', axis = 1)
                    pd.testing.assert_frame_equal(load_section_ID.drop('instantaneous', axis = 1), sections, rtol=1e-5, atol=1e-8)


'''
    def test_import_missing_traj(self): 
        with pytest.raises(FileNotFoundError, match='Trajectory file not found'):
            section_id(topology='data/2f4k.pdb', projection_method='Dihedrals')
    
        def test_import_missing_topo(self): 
        with pytest.raises(FileNotFoundError, match='Topology file not found'):
            section_id(trajectory='data/short_2f4k_f0.xtc', projection_method='Dihedrals')
'''
'''
    def test_import_none(self):
        with pytest.raises(ImportError, match=''):
            section_id(projection_method='Dihedrals')
'''       



class TestSections:
    def test_short_window(self, load_mol):
        with pytest.raises(ValueError, match="`window_size` must be > 1."): 
            section_id(mol=load_mol, projection_method = 'Dihedrals',  window_size=1)

    def test_wrong_method(self, load_mol):
        with pytest.raises(TypeError, match='id_method must be "local" or "global"'):
            section_id(mol=load_mol, projection_method='Dihedrals', id_method='Local')

