from md_intrinsic_dimension import intrinsic_dimension, section_id, secondary_structure_id
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib import cm
from matplotlib.colors import ListedColormap
import matplotlib as mpl
import seaborn as sns
import logging
from moleculekit.molecule import Molecule
import moleculekit.projections.metricrmsd as rmsd
from itertools import product
import mdtraj as md
from matplotlib.ticker import FormatStrFormatter

#build personalised cmap
colors = mpl.colors.ListedColormap(
    name="discrete-bicolor",
    colors=["#e9ff99","#ceff29", "#A5CC21", "#C099F3", "#6100e0", "#4E00B3"])


#set font dimension
plt.rcParams.update({
	'axes.titlesize': 13,
	'axes.labelsize': 13,
	'xtick.labelsize': 11,
	'ytick.labelsize': 11,
	'legend.fontsize': 11,
	'legend.title_fontsize': 13,
	'lines.linewidth' : 1,
	'lines.markersize': 8,
})


folder= 'villin'
topology = '2f4k.pdb'
trajectory = '2f4k'
states = ['u0','u1', 'u2', 'f0', 'f1', 'f2'] 
projection_method='Dihedrals'


##################################### instantaneous ID and RMSD #####################################
ref = Molecule(f'../{folder}/topology')
data = []
for state in states:
    mean_all, mean_last, local_id = intrinsic_dimension(topology=f'../{folder}/{topology}', trajectory=f'../{folder}/{trajectory}_{state}.xtc', id_method = 'local', projection_method=projection_method, verbose=False)
    mol = ref
    mol.read(trajectory+f'_{state}.xtc')
    met=rmsd.MetricRmsd(ref, trajrmsdstr= 'protein and name CA')
    rmsd_values=met.project(mol)
	
    data.append({'trajectory': state,
	  'estimator': 'TwoNN',
	  'mean_all': mean_all,
	  'mean_last': mean_last,
      'local_id': local_id,
	  #'rmsf' : rmsf,
	  'rmsd': rmsd_values})
data = pd.DataFrame(data)
data["folded"] = data["trajectory"].str.startswith("f")

##plots

#merge vertically
fig, ax = plt.subplots(figsize=(10, 15), nrows=3)

# Top: Instantaneous ID vs time
for i, s in enumerate(data['trajectory']):
    local_id = data.loc[i, 'local_id']
    frames = np.arange(len(local_id))
    time_ns = frames / 10.0  # Convert frames to ns
    ax[0].plot(time_ns, local_id, color=colors.colors[i], linewidth=0.8, label=s)

ax[0].set_xlabel("Time (ns)")
ax[0].set_ylabel("ID")
ax[0].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax[0].xaxis.set_major_formatter(FormatStrFormatter('%.1f'))

#Center: ID frequency histogram

all_values = np.concatenate(data['local_id'].values)
bins = np.linspace(all_values.min(), all_values.max(), 30)

for i, s in enumerate(data['trajectory']):
    local_id = data.loc[i, 'local_id']
    ax[1].hist(local_id, bins=bins, alpha=0.7, label=s, 
               color=colors.colors[i], edgecolor='black', linewidth=0.5)

ax[1].set_xlabel("ID")
ax[1].set_ylabel("Number of Frames")
ax[1].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax[1].xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax[1].legend(title="Trajectory")

# Bottom: RMSD vs ID scatter

for i, s in enumerate(states):
    subset = data[data['trajectory'] == s]
    local_id_array = np.concatenate(subset['local_id'].values)
    rmsd_array = np.concatenate(subset['rmsd'].values)
    ax[2].scatter(local_id_array, rmsd_array, 
                  color=colors.colors[i], alpha=0.75, label=s, s=10)

#ax[2].set_box_aspect(1)
ax[2].set_xlabel("ID")
ax[2].set_ylabel("RMSD (Å)")
ax[2].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax[2].xaxis.set_major_formatter(FormatStrFormatter('%.1f'))

for ax_sub, label in zip([ax[0], ax[1], ax[2]], ["A)", "B)", "C)"]):
    ax_sub.annotate(label, xy=(-0.08, 1.03), xycoords='axes fraction',
                    fontsize=16, fontweight='bold', fontfamily='serif',
                    ha='right', va='bottom')

plt.tight_layout()

plt.savefig(f'../../extra/{folder}_combined_instantaneous_RMSD.pdf')
plt.show()


##################################### sections and sec structures #####################################

#sections
all_results = []
for s in states: 
    results = section_id(topology=f'../{folder}/{topology}', trajectory=f'../{folder}/{trajectory}_{s}.xtc', window_size = 15, stride = 3, projection_method=projection_method, id_method='local', verbose=False)
    results["trajectory"] = s  
    all_results.append(results)
results = pd.concat(all_results, ignore_index=True)
results["folded"] = results["trajectory"].str.startswith("f")
results['window'] = results["start"].astype(str) + "-" + results["end"].astype(str)



states = ['u0','u1', 'u2', 'f0', 'f1', 'f2'] 
all_results_ss = []
ss_assignments = []
residue_numbers = None
mol_ref_ss = Molecule(f'/../{folder}/{topology}') #same topology, else change

#sec structures
for s in states:
    results_ss, ss_table = secondary_structure_id(topology=f'../{folder}/{topology}', trajectory=f'../{folder}/{trajectory}_{s}.xtc', mol_ref = mol_ref_ss, simplified=True, projection_method=projection_method, id_method='local' , verbose=False)
    results_ss['trajectory'] = s  
    if residue_numbers is None:
        residue_numbers = ss_table['resid index'].values
    ss_assignments.append(ss_table['sec str type'].values)
    results_ss['trajectory'] = s   
    all_results_ss.append(results_ss)

results_ss = pd.concat(all_results_ss, ignore_index=True)
results_ss["folded"] = results_ss["trajectory"].str.startswith("f")
results_ss['window'] = results_ss["start"].astype(str) + "-" + results_ss["end"].astype(str)

#plot
fig, axes = plt.subplots(ncols=2, figsize=(14,6))  # 2 side-by-side subplots
ax1, ax2 = axes

######################
# Left: sections plot
shift_amount = 0.3
unique_trajectories = results['trajectory'].unique()
n_traj = len(unique_trajectories)
shifts = np.linspace(-shift_amount, shift_amount, n_traj)

for i, trajectory in enumerate(unique_trajectories):
    subset = results[results['trajectory'] == trajectory]
    x_labels = subset['window'].values
    x_positions = np.arange(len(x_labels))
    x_shifted = x_positions + shifts[i]
    y = subset['entire simulation'].values
    y_std = subset['instantaneous'].apply(lambda x: np.std(list(map(float, x.split(',')))))

    ax1.errorbar(x_shifted, y, yerr=y_std, fmt='o',color=colors.colors[i], markeredgecolor='black',ecolor='black', elinewidth=0.8, capsize=3, markersize=8,label=trajectory)


ax1.set_xticks(np.arange(len(x_labels)))
ax1.set_xticklabels(x_labels, rotation=45, fontsize=12)
ax1.set_ylabel("ID", fontsize=13)
ax1.set_xlabel("Window", fontsize=13)
#ax1.set_title("Sections")
ax1.set_ylim(0, max(y)+2)
ax1.set_box_aspect(1)
ax1.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

##plots

# Right: secondary structure plot
shift_amount = 0.3
unique_trajectories = results_ss['trajectory'].unique()
n_traj = len(unique_trajectories)
shifts = np.linspace(-shift_amount, shift_amount, n_traj)

for i, trajectory in enumerate(unique_trajectories):
    subset = results_ss[results_ss['trajectory'] == trajectory]
    x_labels = subset['window'].values
    x_positions = np.arange(len(subset))
    x_shifted = x_positions + shifts[i]
    y = subset['entire simulation'].values
    y_std = subset['instantaneous'].apply(lambda x: np.std(list(map(float, x.split(',')))))

    ax2.errorbar(x_shifted, y, yerr=y_std, fmt='o',color=colors.colors[i], markeredgecolor='black',ecolor='black', elinewidth=0.8, capsize=3, markersize=8,label=trajectory)

all_windows = sorted(subset['window'].unique(), key=lambda w: int(w.split("-")[0]))

xtick_labels = [f"{w} ({s})" for w, s in zip(all_windows, subset['sec str type'])]
ax2.set_xticks(np.arange(len(all_windows)))
ax2.set_xticklabels(xtick_labels, rotation=45, fontsize = 15)
ax2.set_ylabel("ID", fontsize=13)
ax2.set_xlabel("Secondary Structure Element", fontsize=13)
ax2.set_ylim(0, max(y)+2)
ax2.set_box_aspect(1)
ax2.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

for ax, label in zip([ax1, ax2], ["A)", "B)"]):
    ax.annotate(label,xy=(-0.12, 1.02),xycoords='axes fraction',fontsize=16,fontweight='bold',fontfamily='serif', ha='right', va='bottom')

handles, labels = ax2.get_legend_handles_labels()
fig.legend(handles, labels, title="Trajectory", loc='center right', bbox_to_anchor=(1, 0.5))

plt.tight_layout(rect=[0,0,0.92,1]) 
plt.savefig(f'../../extra/{folder}_combined_{projection_method}.pdf', dpi=300)
plt.show()


