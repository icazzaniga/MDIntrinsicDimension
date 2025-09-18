# MDIntrinsicDimension 

MDIntriniscDimension is a Python package for the estimatation of  the intrinsic dimension (ID) of high-dimensional data from, protein molecular dynamics, helping you reduce data to a lower-dimensional manifold while preserving essential information.

It includes three functions and as many analysis modes:    

- `intrinsic_dimension`: ID on whole molecule. 
- `section_id`: ID on fixed sliding window along protein's sequence.    
- `secondary_structure_id`: ID on contiguous secondary structure elements.  

Installation
============

We recommend `uv venv` to create an isolated environment and `uv pip` to install.

To install:

```bash
uv pip install git+https://github.com/icazzaniga/MDIntrinsicDimension.git
```
or:   

```bash
   git clone https://github.com/icazzaniga/MDIntrinsicDimension.git
   uv pip install IntrinsicDimension
```

Quick start
=============
ID can be computed as follows:

```python
    
    #ID of the entire object
    intrinsic_dimension(topology = 'villin/2f4k.pdb', trajectory = 'villin/2f4k_f1.xtc')

    #ID per fixed windows
    section_id(topology = 'villin/2f4k.pdb', trajectory = 'villin/2f4k_f1.xtc')  

    #ID contiguous secondary structure elements
    secondary_structure_id(topology = 'villin/2f4k.pdb', trajectory = 'villin/2f4k_f1.xtc')
```

Any other parameter shared by the functions or specific for each, has a default.  
Please refer to the [documentation]() and the [paper]() for detailed API and tutorials.

Dependencies
------------

These are automatically installed.

- [MoleculeKit](https://software.acellera.com/moleculekit/)
- [scikit-dimension](https://scikit-dimension.readthedocs.io/en/latest/)
- [Numpy](https://numpy.org/)
- [Pandas](https://pandas.pydata.org/)

