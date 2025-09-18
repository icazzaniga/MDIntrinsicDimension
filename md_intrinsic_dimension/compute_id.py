import skdim
import numpy as np
import pandas as pd
from .compute_projections import *
import logging

logger = logging.getLogger(__name__)

def compute_local(projection, estimator = 'TwoNN', last = 100, **id_kwargs):
	'''Computes intrinsic dimension for each frame of the simulation (instantaneous).
	
	Parameters
	----------
	projection : np.ndarray
		Array of projections of shape (frames, features) computed with projections_moleculekit.py
	estimator : str, optional
		Estimator used to compute intrinsic dimension. 
		Options: 'CorrInt', 'DANCo', 'ESS', 'FisherS', 'KNN', 'lPCA', 
		'MADA', 'MiND_ML', 'MLE', 'MOM', 'TLE', 'TwoNN' (default 'TwoNN')
	last : int, optional
		Defines how many frames to consider for mean_last calculation, starting from the end of the simulation (default 100).
    Any other additional keys are passed directly to the chosen estimator’s
    constructor, and should match its parameter names. For example: ``{"estimator": "KNN", "k": 15}``.

	Returns
	-------
	mean_last : float
		Mean of the last `last` local-ID values
	mean_all : float
		Mean of all local-ID values over the trajectory
	local_id : np.ndarray
		Full local-ID time series for each frame, shape (n_frames,)
	'''

	id_estimator = getattr(skdim.id, estimator)(**id_kwargs)#contains only extra parameters
	lid= id_estimator.fit_transform_pw(projection, smooth=True)[1]
	
	mean_last = float(np.mean(lid[-last:]))
	mean_all  = float(np.mean(lid))

	return  mean_all, mean_last, lid



###############################



def compute_global(projection, estimator = 'TwoNN', last = 100, **id_kwargs):
	'''Computes intrinsic dimension for the entire simulation (global).
	
	Parameters
	----------
	projection : np.ndarray
		Array of projections of shape (frames, features) computed with projections_moleculekit.py
	estimator : str, optional
		Estimator used to compute intrinsic dimension. 
		Options: 'CorrInt', 'DANCo', 'ESS', 'FisherS', 'KNN', 'lPCA', 
		'MADA', 'MiND_ML', 'MLE', 'MOM', 'TLE', 'TwoNN' (default 'TwoNN')
	last : int, optional
		Defines how many frames to consider for gid100 calculation, starting from the end of the simulation (default 100)
    Any other additional keys are passed directly to the chosen estimator’s
    constructor, and should match its parameter names. For example: ``{"estimator": "KNN", "k": 15}``.

	Returns
	-------
	gid : float
		Global intrinsic dimension computed over entire trajectory
	gid100 : float
		Global intrinsic dimension computed over last `last` frames

	'''

	id_estimator = getattr(skdim.id, estimator)(**id_kwargs) #contains only extra parameters
	
	gid = id_estimator.fit_transform(projection)
	gid100 = id_estimator.fit_transform(projection[-last:])

	return gid, gid100









