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

def load_molecule(topology = None, trajectory= None, mol = None, verbose=True):
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

