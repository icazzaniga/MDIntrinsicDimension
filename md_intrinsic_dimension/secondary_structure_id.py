from .md_intrinsic_dimension import *
import logging
import numpy as np
import pandas as pd
from moleculekit.molecule import Molecule 
import moleculekit.projections.metricsecondarystructure as mss 
import os 

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
formatter = logging.Formatter('%(name)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.propagate = False

def secondary_structure_id(topology= None, trajectory=None, mol = None, mol_ref=Molecule, simplified=True, projection_method = 'Distances', id_method = 'local', projection_kwargs = None, id_kwargs = None, verbose=True):
    '''
    Computes intrinsic dimension (ID) estimation on contiguous secondary structure elements identified from a protein trajectory.
    This function loads a molecular trajectory, identifies consecutive residues with the same secondary structure assignment (using DSSP via MoleculeKit), 
    and for each secondary structure element computes a molecular projection and estimates the intrinsic dimension using the `scikit-dimension` package.

    Parameters
    ----------
    topology : str, optional
        Path to the topology file (e.g., .pdb, .psf). Required if `mol` is not provided.
    trajectory : str, optional
        Path to the trajectory file (e.g., .dcd, .xtc). Required if `mol` is not provided.
    mol : Molecule, optional
        A pre-loaded MoleculeKit `Molecule` object. If provided, `topology` and `trajectory` are ignored.
    mol_ref : Molecule
        A pre-loaded MoleculeKit  `Molecule` object of ONE FRAME from which DSSP is computed.
    simplified : bool, default=True
        Whether to use the simplified DSSP classification.
        If True (simplified DSSP):
            - 'C' (coil): includes T, S, and loops.
            - 'E' (strand): includes E and B.
            - 'H' (helix): includes H, G, and I.
        If False (full DSSP):
            Residues are assigned full DSSP codes as follows:
                - 'H' : Alpha helix
                - 'B' : Isolated beta-bridge
                - 'E' : Extended strand
                - 'G' : 3₁₀ helix
                - 'I' : π-helix
                - 'T' : Turn
                - 'S' : Bend
                - ' ' (space) : Loops or irregular elements (unassigned)
    projection_method : str or callable, default='Distance'
        Method for generating the molecular projection. Can be one of:
            - 'Distance': pairwise distances or number of contacts between selected atoms.
            - 'Dihedral': backbone or side-chain dihedral angles.
            - MoleculeKit metric class (excluded metrics: Rmsd, SecondaryStructure, TMscore).
    id_method : str, default='local'
        Method for computing intrinsic dimension. One of:
            - 'local': compute frame-wise ID (instantaneous) and average.
            - 'global': compute ID over the entire projection.
    projection_kwargs : dict, optional
        Parameters passed to the projection method. Examples:
            For "Distance":
                - sele : str, atom selection string (default="name CA").
                - step : int, atom subsampling interval (default=1).
                - metric : str, either "distances" or "contacts" (default="distances").
            For "Dihedral":
                - dihedrals : tuple of str, dihedral names (e.g., "phi", "psi", "chi1", ...) (default=("phi", "psi")).
                - sincos : bool, whether to return sin/cos of angles (default=True).
    id_kwargs : dict, optional
        Parameters for intrinsic dimension estimation. Examples:
            - estimator : str, e.g., "TwoNN", "MLE", "KNN", etc. (default="TwoNN").
            - last : int, number of final frames to average over (default=100).
    verbose : bool, default=True
        If True, logs are shown. If False, logs are suppressed.

    Returns
    -------
    results : pandas.DataFrame
        Table with ID results per secondary structure segment.
        Columns: 'window index', 'resids range',sec str type, 'entire simulation', 'last simulation', 'instantaneous'.

    secStr_table : pandas.DataFrame
        Per-residue DSSP assignment.

    Raises
    ------
    FileNotFoundError
        If required topology or trajectory files are missing.
    ValueError
        If the molecule is empty or projection fails.
    TypeError
        If `projection_method` or `id_method` is invalid.
    ImportError
        If a custom projection class cannot be loaded.

    Notes
    -----
    Requires `MoleculeKit` for DSSP projection and `scikit-dimension` for ID estimation.
    '''

        # Configure logger verbosity
    if verbose:
        logger.setLevel(logging.INFO)
    else:
        logger.setLevel(logging.ERROR) #does not show intrinsic_dimension looped
    
    # ----DEFAULT KWARGS PARAMETERS ----
    projection_kwargs = projection_kwargs or {}
    id_kwargs = id_kwargs or {}

    #system preparation and subdivision
    if mol is None:
        if not os.path.isfile(topology):
            raise FileNotFoundError(f'Topology file not found: {topology}')
        if not os.path.isfile(trajectory):
            raise FileNotFoundError(f'Trajectory file not found: {trajectory}')
        
        mol = Molecule(topology, validateElements = False) #ref:PeriodicTable raises error with dummy atoms i. e. M
        mol.read(trajectory)
    
    if mol_ref is None:
        raise FileNotFoundError(f'Missing reference structure for DSSP computation. Please provide a one frame MoleculeKit Molecule object.')
    else:
        if mol_ref.numFrames > 1:
            raise ValueError('ref_mol must be 1 frame long. Please load only the topology in this MoleculeKit Molecule object.')
        if mol.numAtoms != mol_ref.numAtoms:
            raise ValueError(f'mol_ref and mol have a different number of atoms.')
        keys = ['name', 'resid', 'resname', 'chain', 'segid']
        for key in keys:
            if not (mol.get(key) == mol_ref.get(key)).all():
                raise ValueError(f'mol and mol_ref differ in {key}.')

    if simplified == True:
        logger.info(f'Secondary structures considered: Coil (C), Strand (E) and Helix (H).')
    else:
       logger.info("Secondary structure types: Alpha Helix (H), Isolated Beta-Bridge (B), Extended Strand (E), \n "
    "3-Helix (G), 5-Helix (I), Turn (T), Bend (S), Loop/Irregular (' ').")

    logger.info(f'Computing {id_method} Intrinsic Dimension from {projection_method}.')

    met = mss.MetricSecondaryStructure(sel = 'protein', simplified = simplified, integer = False) #integer converts letters to numbers
    
    projection = met.project(mol_ref)

    #projection = met.project(mol) #x: frames; y: ss
    indexes = mol_ref.get('resid', sel='name CA')
    resnames = mol_ref.get('resname', sel='name CA')  #work only on protein, ignore ligands, cofactors, glycans, ...
    

    data = {'resid index': indexes, 'resname': resnames, 'sec str type': projection[0]}
    secStr_table = pd.DataFrame(data)

    start = [] #stores first amino of the sec str
    end = [] #stores last amino of the sec str
    ss_start = secStr_table.iloc[0]['resid index'] #first amino of the protein always start
    ss_type = []
    for i in range(1, len(secStr_table)):
        current = secStr_table.iloc[i]['sec str type']
        previous = secStr_table.iloc[i - 1]['sec str type']
        if current != previous:
            ss_end = secStr_table.iloc[i - 1]['resid index']
            start.append(ss_start)
            end.append(ss_end)
            ss_start = secStr_table.iloc[i]['resid index']
            ss_type.append(previous)
    start.append(ss_start)
    end.append(secStr_table.iloc[-1]['resid index'])
    ss_type.append(secStr_table.iloc[-1]['sec str type'])
    secStr_sequence = list(zip(start, end, ss_type)) 

    #from here compute ID
    results =[]
    for start, end, ss in secStr_sequence:
        if (end - start) < 1: #too short to compute any projection
            logger.warning(f'Skipping segment {start}-{end}: at least two residues per segment are required.')
            continue
        resid_sele = f'resid {start} to {end}'
        window_mol = mol.copy()
        window_mol.filter(resid_sele, _logger=False) ##

        #outname = f'{start:03d}'

        if id_method == 'local':
            all_sim, last, instantaneous = intrinsic_dimension(mol=window_mol, projection_method=projection_method, id_method='local', projection_kwargs=projection_kwargs, id_kwargs=id_kwargs, verbose = False)
        else:  # global
            all_sim, last = intrinsic_dimension(mol=window_mol, projection_method=projection_method, id_method='global', projection_kwargs=projection_kwargs, id_kwargs=id_kwargs, verbose = False)
            instantaneous = []

        results.append({
            'start': start,
            'end': end, 
            'sec str type': ss,
            'window': window_mol.get('resid', 'name CA'),
            'entire simulation': all_sim,
            'last simulation': last, 
            'instantaneous': instantaneous, 
        })        
    results = pd.DataFrame(results)
    return results, secStr_table
