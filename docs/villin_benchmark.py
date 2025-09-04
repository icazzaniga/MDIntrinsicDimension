
#rinominare savefig 


from intrinsic_dimension import intrinsic_dimension, section_id, secondary_structure_id
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap
import logging
from moleculekit.molecule import Molecule
#estimators = ['TwoNN']
#estimators = ['CorrInt', 'ESS', 'FisherS','KNN', 'lPCA', 'MADA', 'MiND_ML', 'MLE', 'MOM', 'TLE', 'TwoNN']
estimators = ['CorrInt', 'DANCo', 'ESS', 'FisherS','KNN', 'lPCA', 'MADA', 'MiND_ML', 'MLE', 'MOM', 'TLE', 'TwoNN']
states = [2, 1]  # trajectory indices
colors = plt.cm.tab10(np.linspace(0, 1, len(states)))  # one color per state

projection_method = 'Dihedrals'  # or 'Distances'
topology = 'examples/villin/2F4K.pdb'
trajectory = f'examples/villin/2F4K_{state}.xtc'

n_estimators = len(estimators)
ncols = 3
fig, axs = plt.subplots(nrows=n_estimators, ncols=ncols, figsize=(15, 5*n_estimators), sharey=False)
axs = np.atleast_2d(axs)

line_handles, line_labels = None, None 

for est_idx, estimator in enumerate(estimators):
    state_ids = []
    means = []
    for state_idx, state in enumerate(states):
        # Compute intrinsic dimension
        mean_all, mean_last, local_id = intrinsic_dimension(topology=topology, trajectory=trajectory, projection_method=projection_method, id_method='local',id_kwargs={'estimator': estimator}, verbose=True)
        state_ids.append(local_id)
        means.append(mean_all)
        frames = np.arange(len(local_id))

        # Column 1: instantaneous ID 
        line, = axs[est_idx, 0].plot(frames, local_id, color=colors[state_idx], lw=0.8, label=state)
        if est_idx == 0:  # grab legend handles only once
            if line_handles is None:
                line_handles, line_labels = [], []
            line_handles.append(line)
            line_labels.append(f"State {state}")

    axs[est_idx, 0].set_ylabel("Local ID")
    if est_idx == n_estimators - 1:
        axs[est_idx, 0].set_xlabel("Frames")

    # Column 2: histogram of instantaneous ID
    all_values = np.concatenate(state_ids)
    bins = np.linspace(all_values.min(), all_values.max(), 30)
    for state_idx, local_id in enumerate(state_ids):
        axs[est_idx, 1].hist(local_id, bins=bins, alpha=0.7, color=colors[state_idx],
                              edgecolor='black', linewidth=0.5)
    axs[est_idx, 1].set_title(estimator)
    axs[est_idx, 1].set_ylabel("Frequency")
    if est_idx == n_estimators - 1:
        axs[est_idx, 1].set_xlabel("Local ID")

    # Column 3: boxplot of local ID per state
bp = axs[est_idx, 2].boxplot(
    state_ids,
    labels=[f"{s}" for s in states],patch_artist=True, showmeans=True,meanprops=dict(marker='*',markerfacecolor='white',markeredgecolor='black', markersize=6)
)

# Color the boxes per trajectory
for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)
    patch.set_edgecolor('black')

for median in bp['medians']:
    median.set_color('black')
axs[est_idx, 2].set_xlabel("Boxplot")
axs[est_idx, 2].set_ylabel("Local ID")  


fig.legend(line_handles, line_labels, title="Trajectory", loc="center left", bbox_to_anchor=(-0.1, 0.5))
plt.tight_layout(rect=[0, 0, 0.95, 0.95])
#plt.savefig('../bin/instantaneous_ID_multi_estimator.png', dpi=300, bbox_inches='tight')
plt.show()
