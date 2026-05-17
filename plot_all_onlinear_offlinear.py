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

fig, axes = plt.subplots(3, 1, figsize=(8, 12))
fig.suptitle('Weighted DFS vs. Haplotype-contiguous DFS', fontsize=14, fontweight='bold')

y_pos = np.arange(len(variant_types))
bar_height = 0.35

def plot_panel_sidebyside(ax, data1, data2, title):
    bars1 = ax.barh(y_pos - bar_height/2, data1, bar_height, color='steelblue', edgecolor='black', linewidth=0.5, label='Weighted DFS')
    bars2 = ax.barh(y_pos + bar_height/2, data2, bar_height, color='lightcoral', edgecolor='black', linewidth=0.5, label='Haplotype-contiguous DFS')
    ax.set_yticks(y_pos)
    ax.set_yticklabels(variant_types)
    ax.invert_yaxis()
    ax.set_title(title, fontsize=11)
    max_val = max(max(data1), max(data2))
    ax.set_xlim(0, max_val * 1.3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.tick_params(left=False, bottom=False)

plot_panel_sidebyside(axes[0], weighted_all, haplo_all, 'All variants')
plot_panel_sidebyside(axes[1], weighted_onlinear, haplo_onlinear, 'On-linear variants')
plot_panel_sidebyside(axes[2], weighted_offlinear, haplo_offlinear, 'Off-linear variants')

# Add panel labels
labels = ['a', 'b', 'c']
for i, ax in enumerate(axes):
    ax.text(-0.15, 1.05, labels[i], transform=ax.transAxes, fontsize=12, fontweight='bold')

from matplotlib.patches import Patch
legend_elements = [Patch(facecolor='steelblue', edgecolor='black', label='Weighted DFS'),
                   Patch(facecolor='lightcoral', edgecolor='black', label='Haplotype-contiguous DFS')]
fig.legend(handles=legend_elements, loc='upper center', ncol=2, bbox_to_anchor=(0.5, 0.02), frameon=False)

plt.tight_layout(rect=(0, 0.04, 1, 1))
plt.savefig('/Users/psalehin/Downloads/weighted_vs_haplo_all_onlinear_offlinear.pdf', bbox_inches='tight')
plt.close()
