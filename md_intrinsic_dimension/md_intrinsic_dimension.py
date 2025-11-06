from .compute_projections import *
from .compute_id import *
import logging
from moleculekit.molecule import Molecule
from moleculekit.projections.projection import Projection
import os 
import importlib 

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
formatter = logging.Formatter('%(name)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.propagate = False



def intrinsic_dimension(topology= None, trajectory=None, mol = None, projection_method = 'Distances', id_method = 'local', projection_kwargs = None, id_kwargs = None, verbose=True):
    '''
    Performs projection of molecular dynamics data followed by intrinsic dimension (ID) estimation.
    This function loads a protein trajectory or a Molecule object from MoleculeKit, computes a projection,
    and estimates the intrinsic dimension using the scikit-dimension package.

    Parameters
    ----------
    topology : str, optional
        Path to the topology file (e.g., .pdb, .psf). Required if `mol` is not provided.
    trajectory : str, optional
        Path to the trajectory file (e.g., .dcd, .xtc). Required if `mol` is not provided.
    mol : Molecule, optional
        A pre-loaded MoleculeKit `Molecule` object. If provided, `topology` and `trajectory` are ignored.
    projection_method : str, callable, or numpy.ndarray, default='Distances'
        Method for generating the molecular projection. Can be one of:
            - 'Distances' : pairwise distances or number of contacts between selected atoms.
            - 'Dihedrals' : selected backbone or side-chain Dihedrals angles.
            - A MoleculeKit metric class name (excluded metrics: Rmsd, SecondaryStructure, TMscore).
            - A MoleculeKit Projection object (projection_kwargs is ignored) 
            - A numpy array, result of a prj.project(mol) operation (trajectory, mol and projection_kwargs are ignored)
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
            - metric : str, either "distances" or "contacts" (default="distances").
            For "Dihedrals"
            - dihedrals : tuple of str, including phi, psi, chi1, .., chi5, omega (default=("psi","phi")).
            - sincos : bool, return sin/cos of angles if True (default=False).
    id_kwargs : dict, optional
        Parameters for intrinsic dimension estimation.
            - estimator : str, name of the estimator from scikit-dimension, including CorrInt, DANCo, ESS, FisherS, KNN, lPCA, MADA, MiND_ML, MLE, MOM, TLE, TwoNN (default="TwoNN").
            - last : int, number of frames to average over starting from the end of the simulation (default=100).
        Additional keys are passed directly to the chosen estimator’s constructor. These should match the estimator’s parameter names in scikit-dimension.
        For example:``{"estimator": "KNN", "k": 15, "last": 200}``
    verbose : bool, default=True
        If True, logging messages are shown. If False, suppress logger output.

    Returns
    -------
    If "local":
	    mean_last : float
		    Mean of the last `last` local-ID values
	    mean_all : float
		    Mean of all local-ID values over the trajectory
	    local_id : np.ndarray
		    Full local-ID time series for each frame, shape (n_frames,)
    
    If "global":
        gid : float
		    Global intrinsic dimension computed over entire trajectory
	    gid100 : float
		    Global intrinsic dimension computed over last `last` frames

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
    Requires `MoleculeKit` for projections and `scikit-dimension` for ID estimation.
    '''
    
    # ----DEFAULT KWARGS PARAMETERS ----
    projection_kwargs = projection_kwargs or {}
    id_kwargs = id_kwargs or {}

    sele = projection_kwargs.get('sele', 'name CA')
    step = projection_kwargs.get('step', 1)
    dihedrals = projection_kwargs.get('dihedrals', ('phi', 'psi'))
    sincos = projection_kwargs.get('sincos', False)
    estimator = id_kwargs.pop('estimator','TwoNN')
    last = id_kwargs.pop('last', int(100))

        # Configure logger verbosity
    if verbose:
        logger.setLevel(logging.INFO)
    else:
        logger.setLevel(logging.CRITICAL + 1)  # effectively disables logger output

    #load Molecule or protein and trajectory
    if mol is None:
        if not os.path.isfile(topology):
            raise FileNotFoundError(f'Topology file not found: {topology}')
        if not os.path.isfile(trajectory):
            raise FileNotFoundError(f'Trajectory file not found: {trajectory}')
        
        mol = Molecule(topology, validateElements = False) #ref:PeriodicTable raises error with dummy atoms i. e. M
        mol.read(trajectory)

    # Determine projection
    builtins = {'Distances': lambda: compute_projections(mol, 'Distances', sele=sele, step=step),
            'Dihedrals': lambda:  compute_projections(mol, 'Dihedrals', dihedrals=dihedrals, sincos=sincos)}
        
    if projection_method in builtins.keys():
        projection = builtins[projection_method]
        projection =projection()
        logger.info(f'Built-in projection "{projection_method}" computed.')
            
    elif isinstance(projection_method, str):
        try:
            module_name = f'moleculekit.projections.metric{projection_method.lower()}'
            class_name = f'Metric{projection_method}'
            module = importlib.import_module(module_name)
            cls = getattr(module, class_name)
            metric = cls(**(projection_kwargs or {}))
            projection = metric.project(mol)
            logger.info(f'Used moleculekit metric projection: {class_name} from {module_name}')
        except Exception as e:
            raise ImportError(
                f'Failed to import or use custom projection class "{class_name}" '
                    f'from module "{module_name}": {e}'
                    )

    elif isinstance(projection_method, Projection):
        projection = projection_method.project(mol)
        logger.info("Projecting using the provided Projection object")

    elif isinstance(projection_method, np.ndarray):
        projection = projection_method
        logger.info("Using the results of a project() operation provided as an array of shape "+str(projection.shape))

    else:
        raise TypeError('projection_method must be a string referring to a MoleculeKit Projection class, an array, or a Projection object. If string, the first letter of each word should be capitalized.')

#AGGIUNGERE CHECK ARRAY

    # ID estimation mapping
    if id_method == 'local':
        logger.info(f'Computing {id_method} intrinsic dimension using estimator "{estimator}" (last simulation section = {last} frames).')
        out = compute_local(projection=projection, estimator=estimator, last=last, **id_kwargs)
    if id_method == 'global':
        logger.info(f'Computing {id_method} intrinsic dimension using estimator "{estimator}" (last simulation section = {last} frames).')
        out = compute_global(projection=projection, estimator=estimator, last=last, **id_kwargs)
    else:
        raise TypeError(
            f'id_method must be "local" or "global", got {id_method} instead.'
        )

    return out
