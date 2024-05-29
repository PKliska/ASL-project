import sys
from pathlib import Path
import math

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

import pandas as pd

measurements_file = Path(sys.argv[1]).resolve(strict=True)
plot_type = sys.argv[2]

df = pd.read_csv(measurements_file)

if has_diff_matrix_dims := len(df['Matrix dimension'].unique()) > 1:
    raise Exception("Error: Trying to plot single heatmap for diff matrix sizes (they need to be the same)")

matrix_dimension = df['Matrix dimension'][0]

df['Performance'] = df['FLOPs'] / df['Cycles']
df['Runtime'] =  df['Cycles'] / (3.4*10**9)

data_perf = df.pivot_table('Performance', 'Block size', 'Timestamps', fill_value=float("nan"))
data_runtime = df.pivot_table('Runtime', 'Block size', 'Timestamps', fill_value=float("nan"))

# Dimensions for the x and y axes
x_labels = df['Timestamps'].unique()
y_labels = df['Block size'].unique()

# Create the heatmap
plt.figure(figsize=(25, 15))
if plot_type.startswith("perf"):
    sns.heatmap(data_perf, annot=True, fmt=".2f", cmap="YlOrRd", xticklabels=x_labels, yticklabels=y_labels)
else:
    # reverse colors cuz lower (runtime) is better
    sns.heatmap(data_runtime, annot=True, fmt=".2f", cmap="YlOrRd_r", xticklabels=x_labels, yticklabels=y_labels)

# Add labels and title
plt.xlabel("Timestamp size")
plt.ylabel("Block size")
if plot_type.startswith("perf"):
    plt.title(f"Ya moma's performance (flop/cycle) for matrix_dim={matrix_dimension}")
else:
    plt.title(f"Ya moma's runtime (sec) for matrix_dim={matrix_dimension}")

# Show the heatmap
plt.show()

base_dir = measurements_file.parent
if plot_type.startswith("perf"):
    heatmap_path = base_dir / 'heatmap_perf.png'
else:
    heatmap_path = base_dir / 'heatmap_runtime.png'

plt.savefig(heatmap_path, dpi=300)
print(f"Plot saved as '{heatmap_path}'")
