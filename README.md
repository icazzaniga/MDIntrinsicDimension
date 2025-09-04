# Intrinsic Dimension Computation over MD Trajectories.

This package estimates the intrinsic dimension (ID) of high-dimensional data from complex biological systems (e.g., protein molecular dynamics), helping you reduce data to a lower-dimensional manifold while preserving essential information.

It includes three primary functions:
- `intrinsic_dimension`: Computes the ID over the entire MD system.
- `section_id`: Computes ID in sliding, fixed-length windows along the protein chain.
- `secondary_structure_id`: Compute ID separately for each secondary structure element.

## Getting started

This package uses pip-uv to create an isolated environment. <br>
Installation:
```sh
uv pip install git+git@github.com:icazzaniga/IntrinsicDimension.git
```
To activate the environment, from within the main IntrinsicDimension directory:
```sh
source env/bin/activate
```
## Dependencies

This package relies on:
- [MoleculeKit projections](https://software.acellera.com/moleculekit/projections.html) for the computation of projections from MD trajectories. 

- [scikit-dimension](https://scikit-dimension.readthedocs.io/en/latest/) to select the ID estimators.

## Warnings and notes

> **Warning:** <br> The ID matrix is of form `n_frames x m_features`, with `n > 100` and `m > 1`. <br> Only the following MoleculeKit projection classes are supported:
> - "Coordinate"         
> - "Dihedrals"        
> - "Distance"
> - "Plumed2"            
> - "Sasa"    
> - "Shell"     
> - "SphericalCoordinate".

> **Warning:** <br> `section_id` and `secondary_structure_id` require respectively `window_size` and the secondary structure element to be longer than 1 residue to be able to compute projections and consequently ID. In case of secondary structures <= 1 ID is not computed but the structure is included in the DataFrames. 

> **Note on Dihedrals and Distances** <br> These two projections are not directly called from MoleculeKit; instead, this package uses custom functions that accept additional parameters for more flexible analysis.

> **Note on ID estimators:** <br>While all the scikit-dimension available estimators can, in principle, be used in this package, the complexity associated to a MD simulation can lead some of the estimators to failure. <br> We set `TwoNN` estimator as default as it has proven to be one of the most robust.

> **Potential file format issues:** <br> MoleculeKit supports many trajectory and topology formats but may raise some parsing problems. <br> We suggest to work with `PDB` for topologies and `dcd` or `xtc` for trajectories. 

## Usage
A [test notebook](examples/test.ipynb) is available in the `examples/` directory, including:
- basic functions usage and data handling
- topology and trajectory loading
- **kwargs usage and examples
- plotting.

