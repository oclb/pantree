import matplotlib.pyplot as plt
import numpy as np

# Data calculated from the figure
# All = Small + Large
# On-linear = All - Off-linear

variant_types = ['SNPs', 'MNPs', 'Insertions', 'Deletions', 'Replacements', 'Duplications', 'Inversions']

# Weighted DFS (blue)
weighted_all = [21682923, 660525, 3245665, 3631665, 336546, 50, 860]
weighted_offlinear = [1646297, 216065, 756002, 722086, 131769, 1, 192]
weighted_onlinear = [20036626, 444460, 2489663, 2909579, 204777, 49, 668]

# Haplotype-contiguous DFS (red)
haplo_all = [21614394, 710924, 3093687, 3753048, 388830, 51, 865]
haplo_offlinear = [1577758, 265275, 604024, 843469, 183628, 0, 350]
haplo_onlinear = [20036636, 445649, 2489663, 2909579, 205202, 51, 515]

fig, axes = plt.subplots(3, 2, figsize=(10, 10))
fig.suptitle('Weighted DFS vs. haplotype-contiguous DFS', fontsize=14, fontweight='bold')

y_pos = np.arange(len(variant_types))

def plot_panel(ax, data, title, color):
    bars = ax.barh(y_pos, data, color=color, edgecolor='black', linewidth=0.5)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(variant_types)
    ax.invert_yaxis()
    ax.set_title(title, fontsize=11)
    for i, v in enumerate(data):
        ax.text(v + max(data)*0.01, i, str(v), va='center', fontsize=8)
    ax.set_xlim(0, max(data) * 1.25)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.tick_params(left=False, bottom=False)

# Row 1: All variants
plot_panel(axes[0, 0], weighted_all, 'All variants', 'steelblue')
plot_panel(axes[0, 1], haplo_all, 'All variants', 'lightcoral')

# Row 2: On-linear variants
plot_panel(axes[1, 0], weighted_onlinear, 'On-linear variants', 'steelblue')
plot_panel(axes[1, 1], haplo_onlinear, 'On-linear variants', 'lightcoral')

# Row 3: Off-linear variants (Non-GRCh38)
plot_panel(axes[2, 0], weighted_offlinear, 'Off-linear variants', 'steelblue')
plot_panel(axes[2, 1], haplo_offlinear, 'Off-linear variants', 'lightcoral')

# Add panel labels
labels = ['a', 'b', 'c', 'd', 'e', 'f']
for i, ax in enumerate(axes.flat):
    ax.text(-0.15, 1.05, labels[i], transform=ax.transAxes, fontsize=12, fontweight='bold')

from matplotlib.patches import Patch
legend_elements = [Patch(facecolor='steelblue', edgecolor='black', label='Weighted DFS'),
                   Patch(facecolor='lightcoral', edgecolor='black', label='Haplotype-contiguous DFS')]
fig.legend(handles=legend_elements, loc='upper center', ncol=2, bbox_to_anchor=(0.5, 0.02), frameon=False)

plt.tight_layout(rect=(0, 0.04, 1, 1))
plt.savefig('/Users/psalehin/Downloads/weighted_vs_haplo_all_onlinear_offlinear.pdf', bbox_inches='tight')
plt.close()
