# MDIntrinsicDimension 

MDIntriniscDimension is a Python package for the estimatation of  the intrinsic dimension (ID) of high-dimensional data from, protein molecular dynamics, helping you reduce data to a lower-dimensional manifold while preserving essential information.

It includes three functions and as many analysis modes:    

- `intrinsic_dimension`: ID on whole molecule. 
- `section_id`: ID on fixed sliding window along protein's sequence.    
- `secondary_structure_id`: ID on contiguous secondary structure elements.  

Installation
------------

We recommend `uv venv` to create an isolated environment and `uv pip` or `uv sync`  to install:

```bash
uv pip install git+https://github.com/icazzaniga/MDIntrinsicDimension.git
```
or:   

```bash
git clone https://github.com/icazzaniga/MDIntrinsicDimension.git
uv sync
```

Documentation
------------

It is available online at https://giorginolab.github.io/MDIntrinsicDimension/ .


Quick start
------------
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
Please refer to the [documentation](https://giorginolab.github.io/MDIntrinsicDimension/) and the [preprint](http://arxiv.org/abs/2511.13550) for detailed API and tutorials.



Example notebooks
-----------------


The documentation includes example notebooks that demonstrate the package’s features and reproduce (and extend) the figures from the paper. To run these notebooks, you have two options:

1.	Use the original DESRES MD trajectories of villin and NTL9 by [Lindorff-Larsen et al.](http://science.sciencemag.org/content/334/6055/517). These must be downloaded, pre-processed, and placed in the `examples` directory.
2.	Use the provided pre-computed data tables, available in CSV format, which require no trajectory downloads.

The pre-processing steps are described at the top of the [Villin Headpiece](docs/villin_plots.ipynb) and [NTL9](docs/NTL9_plots.ipynb) notebooks. The same notebooks also appear in the online [documentation](https://giorginolab.github.io/MDIntrinsicDimension/) under the Examples section, with some code hidden for readability. The full source code is included in the notebooks themselves and associated files.



Dependencies
------------

They are automatically installed.

- [MoleculeKit](https://software.acellera.com/moleculekit/)
- [scikit-dimension](https://scikit-dimension.readthedocs.io/en/latest/)



Citation
------------

Consider citing the package through the corresponding paper (under review):

> I. Cazzaniga, T. Giorgino. MDIntrinsicDimension: Dimensionality-Based Analysis of Collective Modes of Macromolecules from Molecular Dynamics Trajectories. 


Notes for developers
------------

To install along with other dependencies used for development (e.g. Jupyter), use

    uv sync --all-groups


To build the docs, use

    uv run --only-group docs make -C docs html
    open docs/_build/html/index.html

