Quick start
=============
Basic Usange
----------------

First, activate the environment:

.. code-block:: bash

   source env/bin/activate

Then, in Python, import the package:

.. code-block:: python

    from intrinsic_dimension import intrinsic_dimension, section_id, secondary_structure_id

At this point, ID can be computed as follows:

.. code-block:: python
    
    #ID of the entire object
    intrinsic_dimension(topology = 'villin/2F4K.pdb', trajectory = 'villin/2F4K_F.xtc')

    #ID of the object divided per windows
    section_id(topology = 'villin/2F4K.pdb', trajectory = 'villin/2F4K_F.xtc')  

    #ID of the object divided bu secondary structure
    secondary_structure_id(topology = 'villin/2F4K.pdb', trajectory = 'villin/2F4K_F.xtc')

As any other parameter shared by the three functions or specific for each, has a default.
These parameters are: 

* :mark:`projection_method`, default "Distances" 
* :mark:`id_method`, default "local" 
* :mark:`projection_kwargs`, including:      
   
    * ``sele``, default "name CA"    
    * ``step``, default 1     
    * ``metric``, default "distances" 

* :mark:`id_kwargs`, extra parameters including:  
    
    * ``estimator``, allows to chose the estimator of desire (default "TwoNN").  
    * ``last``, for more precise results, all the functions in the package allow the computation of ID on the last part of the trajectory (default "100").   

And, in case of **section_id**: 

* :mark:`window_size`, default 10
* :mark:`stride`, default 1

Wheras, for **secondary_structure_id**:

* :mark:`simplified`, default True
  
File Format Compatibility 
-------------------------------------
The package supports many file formats, but we recommend using:

- `pdb` for topology
- `xtc` or `dcd` for trajectory

.. attention::

    **scikit-dimension** requires **at least 101 frames** of trajectory to repourpose a global estimator as local as the default neighbourhood is of 100 elements.


Non-basic usage
----------------
It is possible to use different parameters than the default ones defined above, for example:

Load the moleculekit `Molecule` object outside the function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. code-block:: python
    
    #create the Molecule object by loading topology and trajectory
    mol = Molecule('villin/2F4K.pdb')
    mol.read('villin/2F4K_F.xtc')

    #call the ID function
    intrinsic_dimension(mol = mol) #topology and trajectory ignored if mol is present

It is also possible to change the default parameters as follows:

Change :mark:`projection_method` and :mark:`projection_kwargs`.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Projections are used as preliminary step for the reduction of the total dimension number. 
Several projection types can be computed, relying on Moleculekit's projections avaialbility.
The selected projection must be called as a string with the first letter in upper case:

.. code-block:: python

    intrinsic_dimension(topology = 'villin/2F4K.pdb', trajectory = 'villin/2F4K_F.xtc', 
            projection_method = 'Dihedral')   

:mark:`projection_kwargs` is a optional input dictionary containing of parameters required by the selected :mark:`projection_method`.    
Default values are provided for methods like ``Distances`` and ``Dihedrals`` (see *Important* below), which are not handled directly via MoleculeKit.  
In particular:
   
    * ``sele`` , default "name CA"    
    * ``step`` , default 1     
    * ``metric`` , default "distances" 
    * ``dihedrals`` , default ("psi", "phi")
    * ``sincos`` , default False

.. code-block:: python

    #define new parameters for projection method "Coordianate"
    proj = {'atomsel':'name CA','refmol':mol} 
    intrinsic_dimension(topology = 'villin/2F4K.pdb', trajectory = 'villin/2F4K_F.xtc', 
    projection_method = 'Coordinate', projection_kwargs=proj)   

.. attention::
    The :mark:`projection_method` string must have the first letter of each word in upper case, the remaining in lower case.
    

.. important:: 
    The ID matrix must be of shape **n_frames x m_features** with n > 100 and m > 1.
    Accordingly, only the following MoleculeKit projection classes are supported:

    * "Coordinate",  
    * "Dihedral",
    * "Distance", 
    * "Fluctuation",
    * "Gyration",
    * "Plumed2",  
    * "Sasa",  
    * "Shell",  
    * "SphericalCoordinate".

    ``Distances`` and ``Dihedrals`` (plural) are modified functions of the MoleculeKit package that accept additional parameters for a more flexible analysis.
    The singular form (``Distance`` and ``Dihedral``), still allow to use the original projection. 

Change :mark:`id_method` and :mark:`id_kwargs`.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
:mark:`id_method` includes ``local``, ``global`` (default "local"). 

In case of ``local`` ID estimation, the estimator identifies sub-regions of the dataset based on a shared local feature, i.e. time, on which ID is computed.     
``global`` ID estimation consists in the computation of a single value of ID for the entire system, it might not be the best choice in case of highily complex systems. 
 
.. code-block:: python
    
    #global or local ID method
    intrinsic_dimension(topology = 'villin/2F4K.pdb', trajectory = 'villin/2F4K_F.xtc', 
            id_method = 'global')   
  

:mark:`id_kwargs` is an optional input dictionary containing:
    
    * ``estimator`` , allows to chose the estimator of desire,default "TwoNN" (see *Important* below).  
    * ``last`` , for more precise results, all the functions in the package allow to exclude the initial part of the trajectory as the molecule might not be stable, slicing from the end of the simulation (default "100").   

.. code-block:: python
    
    #define id_kwargs parameters 
    estimator = {'estimator': 'lPCA', 'last': '100'} 
    intrinsic_dimension(topology = 'villin/2F4K.pdb', trajectory = 'villin/2F4K_F.xtc',
             id_kwargs= estimator)   


.. important:: 
    While all the scikit-dimension available `estimators <https://scikit-dimension.readthedocs.io/en/latest/api.html>`_ can, in principle, be used in this package, the complexity associated to a MD simulation can lead some of the estimators to failure.
    We suggest using **TwoNN** estimator (default) as it has proven to be one of the most robust.












