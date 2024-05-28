import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

import pandas as pd

df = pd.read_csv('performance_metrics.csv')

if has_diff_matrix_dims := len(df['Matrix dimension'].unique()) > 1:
    raise Exception("Error: Trying to plot single heatmap for diff matrix sizes (they need to be the same)")

matrix_dimension = df['Matrix dimension'][0]

df['Performance'] = df['FLOPs'] / df['Cycles']
df['Runtime'] =  df['Cycles'] /(3.4*10**9)

data_perf = df.pivot_table('Performance', 'Block size', 'Timestamps', fill_value=0)
data_runtime = df.pivot_table('Runtime', 'Block size', 'Timestamps', fill_value=0)


# Dimensions for the x and y axes
x_labels = [1,2,5,10,25]
y_labels = [27,32,36,48,54,72,81,96, 108]

# Create the heatmap
plt.figure(figsize=(8, 6))
if sys.argv[1].startswith("perf"):
    sns.heatmap(data_perf, annot=True, fmt=".1f", cmap="YlOrRd", xticklabels=x_labels, yticklabels=y_labels)
else:
    # reverse colors cuz lower (runtime) is better
    sns.heatmap(data_runtime, annot=True, fmt=".1f", cmap="YlOrRd_r", xticklabels=x_labels, yticklabels=y_labels)

# Add labels and title
plt.xlabel("Timestamp size")
plt.ylabel("Block size")
if sys.argv[1].startswith("perf"):
    plt.title(f"Ya moma's performance (flop/cycle) for matrix_dim={matrix_dimension}")
else:
    plt.title(f"Ya moma's runtime (sec) for matrix_dim={matrix_dimension}")


# Show the heatmap
plt.show()
if sys.argv[1].startswith("perf"):
    plt.savefig("heatmap_perf.png", dpi=300)
    print(f"Plot saved as '{Path('./heatmap_perf.png').resolve(strict=True)}'")
else:
    plt.savefig("heatmap_runtime.png", dpi=300)
    print(f"Plot saved as '{Path('./heatmap_runtime.png').resolve(strict=True)}'")
