import moleculekit.projections.metricdistance as distance
from moleculekit.projections.metricdihedral import MetricDihedral, Dihedral
import numpy as np


def compute_projections(mol, projection_method, **kwargs):
    '''
    Computes molecular projections (features) from a MoleculeKit `Molecule` object.

    Supports distance-based and dihedral-based features for downstream analysis such as 
    dimensionality reduction or intrinsic dimension estimation.

    Parameters
    ----------
    mol : moleculekit.molecule.Molecule
        MoleculeKit object containing atomic structure and trajectory.
    projection_method : str
        Type of projection to compute:
            - 'Distances' : pairwise distances between selected atoms.
            - 'Dihedrals' : specified backbone or side-chain dihedral angles.
    **kwargs : dict
        Extra arguments specific to the projection method:

        For 'Distances':
            sele : str, default='protein and name CA'
                Atom selection string (VMD format).
            step : int, default=1
                Subsampling step over selected atoms.
            metric : str, default='distances'
                Either 'distances' or 'contacts'.

        For 'Dihedrals':
            dihedrals : tuple of str, default=('phi', 'psi')
                Dihedral angles to compute.
            sincos : bool, default=False
                If True, return sine and cosine of angles instead of degrees.

    Returns
    -------
    projection : np.ndarray
        2D array of shape (n_frames, m_features) with computed features per frame.
        NOTE: n > 100, m > 1

    Raises
    ------
    ValueError
        If input molecule is empty or projection method is invalid.
    '''
    #check if non empty
    num_frames = getattr(mol, 'numFrames', None)
    num_atoms = getattr(mol, 'numAtoms', None)
    if num_frames is None or num_frames == 0:
        raise ValueError('Provided Molecule contains no trajectory frames.')
    if num_atoms is None or num_atoms == 0:
        raise ValueError('Provided Molecule contains no atoms.')
    

    if projection_method == 'Distances':
        sele = kwargs.get('sele', 'protein and name CA')   ###########
        step = kwargs.get('step', 1)
        metric_type = kwargs.get('metric', 'distances') #default is distances
        if metric_type not in ('distances', 'contacts'):
            raise ValueError(f'Invalid metric type: {metric_type}. Use 'distances' or 'contacts'.')
        all_atoms = mol.atomselect(sele, indexes=True)
        
        atoms = all_atoms[0::step]
        dim = len(atoms) * (len(atoms) - 1) // 2
        pairs = np.zeros((dim, 2), dtype=int)
        k = 0
        for i in range(len(atoms)):
            for j in range(i + 1, len(atoms)):
                pairs[k] = [atoms[i], atoms[j]]
                k += 1
        
        sel1 = [f'index {i}' for i in pairs[:, 0]]
        sel2 = [f'index {j}' for j in pairs[:, 1]]

        
        #sel1 = 'protein and name CA'
        #sel2 = 'protein and name CA'
        met = distance.MetricDistance(sel1=sel1, sel2=sel2,
                                    metric=metric_type, periodic='selections') #also contacts
        projection = met.project(mol)
        return projection

    elif projection_method == 'Dihedrals':
        dihedrals = kwargs.get('dihedrals', ('phi', 'psi'))
        sincos = kwargs.get('sincos', False)
        angles = Dihedral.proteinDihedrals(mol=mol, dih=dihedrals)
        met = MetricDihedral(dih=angles, sincos=sincos)
        projection = met.project(mol)
        return projection
        



#chi1 ARG ASN ASP CYS GLN GLU HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL
#chi2 ARG ASN ASP GLN GLU HIS ILE LEU MET PHE PRO TRP TYR
#chi3 ARG GLN GLU LYS MET
#chi4 ARG LYS
#chi5 ARG






