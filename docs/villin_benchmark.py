



from intrinsic_dimension import intrinsic_dimension, section_id, secondary_structure_id
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap
import logging
from moleculekit.molecule import Molecule

states = [1]
colors = plt.cm.jet(np.linspace(0, 1, len(states)))
estimators = ['CorrInt', 'ESS', 'FisherS','KNN', 'lPCA', 'MADA', 'MiND_ML', 'MLE', 'MOM', 'TLE', 'TwoNN']

mol_proj = 'Distance'
proj_kwargs={'skip': 3}

mol = Molecule('examples/villin/2F4K.pdb')
mol.read('examples/villin/2F4K_1.xtc')

# Grid size (make a tall column with 2 per row)
ncols = 2
nrows = int(np.ceil(len(estimators) / ncols))

fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(6 * ncols, 3 * nrows), sharey=True)
axs = np.atleast_1d(axs).flatten()  # Flatten for easy iteration

for idx, estimator in enumerate(estimators):
    ax = axs[idx]
    mean_all, mean_last, local_id = intrinsic_dimension(
        mol=mol,
        projection_method=mol_proj,
        id_method='local',
        id_kwargs={'estimator': estimator},
        projection_kwargs = proj_kwargs,
        verbose=True
    )
    frames = np.arange(len(local_id))
    ax.plot(frames, local_id, color=colors[0], linewidth=0.8, label=f"State 1")

    ax.set_title(f"{estimator}")
    ax.set_xlabel("Frames")
    if idx % ncols == 0:  # only left column gets ylabel
        ax.set_ylabel("Local ID")

# Remove unused subplots if any
for j in range(len(estimators), len(axs)):
    fig.delaxes(axs[j])

# Optional legend outside
# handles, labels = ax.get_legend_handles_labels()
# fig.legend(handles, labels, bbox_to_anchor=(1.05, 0.5), loc='center left', title='State')

plt.tight_layout(rect=[0, 0, 0.9, 0.95])
plt.savefig("villin_benchmark_2.png", dpi=300, bbox_inches='tight')
plt.show()

