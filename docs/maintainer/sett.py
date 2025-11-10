from md_intrinsic_dimension import intrinsic_dimension, section_id, secondary_structure_id
import numpy as np
import pandas as pd
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap
import matplotlib as mpl
import seaborn as sns
import logging
from moleculekit.molecule import Molecule 
import moleculekit.projections.metricrmsd as metricrmsd
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
from itertools import product
from matplotlib.lines import Line2D
from pathlib import Path


#build personalised cmap
colors = mpl.colors.ListedColormap(
    name="discrete-bicolor",
    colors=["#e9ff99","#ceff29", "#A5CC21", "#C099F3", "#6100e0", "#4E00B3"])
cmap = plt.get_cmap('jet')

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

print('Settings loaded.')