from .intrinsic_dimension import *
import logging
import numpy as np
from moleculekit.molecule import Molecule 
import os 
import importlib 

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
formatter = logging.Formatter('%(name)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.propagate = False


def section_intrinsic_dimension(topology=None, trajectory=None, mol=None, window_size=10, stride=1, projection_method='Distances', id_method='local', projection_kwargs=None, id_kwargs=None, verbose = True):
    '''
        Performs projection of molecular dynamics data followed by intrinsic dimension (ID) estimation.
    This function loads a protein trajectory in MD, computes a projection (e.g., distances or dihedrals),
    and estimates the intrinsic dimension using the `scikit-dimension` package.

    Parameters
    ----------
    topology : str, optional
        Path to the topology file (e.g., .pdb, .psf). Required if `mol` is not provided.
    trajectory : str, optional
        Path to the trajectory file (e.g., .dcd, .xtc). Required if `mol` is not provided.
    mol : Molecule, optional
        A pre-loaded MoleculeKit `Molecule` object. If provided, `topology` and `trajectory` are ignored.
    window_size : int, default=10, number of amino acids to be considered per segment.
    stride : int, default=1, number of amino acids between one window and the following.
    projection_method : str or callable, default='Distances'
        Method for generating the molecular projection. Can be one of:
            - 'Distances' : pairwise distances or number of contacts between selected atoms.
            - 'Dihedrals' : selected backbone or side-chain dihedral angles.
            - Any custom MoleculeKit metric projection class name (e.g., 'Rmsd', 'Sasa').
    id_method : str, default='local'
        Method for computing intrinsic dimension. One of:
            - 'local' : compute frame-wise ID (instantaneous) and averaged.
            - 'global' : compute ID over the entire projection.
    projection_kwargs : dict, optional
        Parameters passed according to the projection method. 
        Defaults:
            For "Distances"
            - sele : str, atom selection string (default="name CA").
            - step : int, subsampling interval (default=1).
            - metric : str, either "distances" or "contacts" (default="distances");
            For "Dihedrals"
            - dihedrals : tuple of str, including phi, psi, chi1, .., chi5, omega (default=("psi","phi")).
            - sincos : bool, return sin/cos of angles if True (default=True).
    id_kwargs : dict, optional
        Parameters for intrinsic dimension estimation. Examples:
            - estimator : str, name of the estimator from scikit-dimension, including CorrInt, DANCo, ESS, FisherS, KNN, lPCA, MADA, MiND_ML, MLE, MOM, TLE, TwoNN (default="TwoNN")
            - last : int, number of frames to average over starting from the end of the simulation (default=100).
    verbose : bool, default=True
        If True, logging messages are shown. If False, suppress logger output.

    Returns
    -------
    results : DataFrame 
        columns include "window:index", resids", "entire simulation", "last simulation" and "instantaneous". 
    
    Raises
    ------
    FileNotFoundError
        If required topology or trajectory files are missing.
    ValueError
        If molecule is empty or projection fails.
    TypeError
        If `projection_method` or `id_method` is invalid.
    ImportError
        If custom projection class cannot be loaded.

    Notes
    -----
    Requires `moleculekit` for projections and `scikit-dimension` for ID estimation.
    '''
    if mol is None:
        if not os.path.isfile(topology):
            raise FileNotFoundError(f'PDB file not found: {topology}')
        if not os.path.isfile(trajectory):
            raise FileNotFoundError(f'DCD file not found: {trajectory}')
        
        mol = Molecule(topology, validateElements = False) #ref:PeriodicTable raises error with dummy atoms i. e. M
        mol.read(trajectory)

    # Configure logger verbosity
    if verbose:
        logger.setLevel(logging.INFO)
    else:
        logger.setLevel(logging.CRITICAL + 1)  # effectively disables logger output

    resids = mol.get('resid', sel='protein')  #work only on protein, ignore ligands, cofactors, glycans, ...
    resids = np.unique(resids) #one number per resid instead of per atom

    total_resids = len(resids)
    windows_number = (total_resids - window_size) // stride
    extra_aa = (total_resids - window_size) % stride
    if extra_aa != 0:
        logger.info(f'Protein has {total_resids} amino acids. Slicing in {windows_number} windows of {window_size} amino acids each and {stride} aminos stride.')
        logger.info(f'Last {extra_aa} amino acids will be ingored.')
    logger.info(f'Computing {id_method} Intrinsic Dimension from {projection_method}.')
    results=[]
    for i in range(0, len(resids) - window_size + 1, stride):
        resid_window = resids[i:i + window_size]
        resid_sele = f'protein and resid {' '.join(map(str, resid_window))}'
        window_mol = mol.copy()
        window_mol.filter(resid_sele, _logger=False)

        outname = f'{i+1:03d}'
        if id_method == 'local':
            all_sim, last, instantaneous = intrinsic_dimension(mol=window_mol, projection_method=projection_method, id_method='local', projection_kwargs=projection_kwargs, id_kwargs=id_kwargs, verbose = False)
        else:  # global
            all_sim, last = intrinsic_dimension(mol=window_mol, projection_method=projection_method, id_method='global', projection_kwargs=projection_kwargs, id_kwargs=id_kwargs, verbose = False)
            instantaneous = []

        results.append({
            'window_index': outname,
            'resids': (min(resid_window), max(resid_window)),
            'entire simulation': all_sim,
            'last simulation': last,
            'instantaneous': instantaneous,
        })
    results = pd.DataFrame(results)

    return results


