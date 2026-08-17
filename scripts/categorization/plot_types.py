import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import configparser

# Load configuration
config = configparser.ConfigParser()
config.read('config.ini')

# Get paths from config
base_path = config['categorization']['base_path']
output_dir = os.path.join(base_path, config['categorization']['output_dir'])

# Colours per measurement type, in column order. ColorBrewer PuOr, which is
# colourblind-safe; listed explicitly so the assignment cannot drift.
colors = [c.strip() for c in config['categorization']['colors'].split(',')]

# Set font
plt.rcParams['font.family'] = 'DejaVu Sans'

# Read data
df = pd.read_csv(os.path.join(output_dir, 'types.txt'), sep='\t', index_col=0)

measurement_types = df.columns.tolist()
color_map = {mtype: colors[i % len(colors)] for i, mtype in enumerate(measurement_types)}

# Version 1: Horizontal bars
fig, ax = plt.subplots(figsize=(12, 6))
left = np.zeros(len(df))
for mtype in measurement_types:
    values = df[mtype].fillna(0).values
    bars = ax.barh(df.index, values, left=left, color=color_map[mtype],
                    edgecolor='white', linewidth=0.5, label=mtype)
    for i, (bar, val) in enumerate(zip(bars, values)):
        if val > 0:
            x_pos = left[i] + val / 2
            ax.text(x_pos, i, f'{int(val)}', ha='center', va='center',
                   fontsize=9, fontweight='bold', color='black')
    left += values

ax.set_xlabel('')
ax.set_ylabel('')
ax.set_xticks([])
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.tick_params(axis='y', labelsize=13)
ax.legend(loc='upper right', frameon=False, fontsize=13)
plt.tight_layout()
plt.savefig(os.path.join(output_dir, 'types_barplot_horizontal.png'),
            dpi=300, bbox_inches='tight')
plt.close()

# Version 2: Vertical bars
fig, ax = plt.subplots(figsize=(10, 7))
bottom = np.zeros(len(df))
for mtype in measurement_types:
    values = df[mtype].fillna(0).values
    bars = ax.bar(df.index, values, bottom=bottom, color=color_map[mtype],
                   edgecolor='white', linewidth=0.5, label=mtype)
    for i, (bar, val) in enumerate(zip(bars, values)):
        if val > 0:
            y_pos = bottom[i] + val / 2
            ax.text(i, y_pos, f'{int(val)}', ha='center', va='center',
                   fontsize=9, fontweight='bold', color='black')
    bottom += values

ax.set_xlabel('')
ax.set_ylabel('')
ax.set_yticks([])
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.tick_params(axis='x', labelsize=13)
ax.legend(loc='upper right', frameon=False, fontsize=13)
plt.tight_layout()
plt.savefig(os.path.join(output_dir, 'types_barplot_vertical.png'),
            dpi=300, bbox_inches='tight')
plt.close()

print("Barplots saved to", output_dir)
