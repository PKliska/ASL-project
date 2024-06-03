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
    raise Exception("Error: Trying to plot single heatmap for diff matrix sizes (they need to be the same)!")
matrix_dimension = df['Matrix dimension'][0]

# if has_diff_timestamp_dims := len(df['Timestamp'].unique()) > 1:
#     raise Exception("Error: Trying to plot single heatmap for diff timestamps (they need to be the same)!")
# timestamp = df['Timestamp'][0]
timestamp = "max"

df['Performance'] = df['FLOPs'] / df['Cycles']
df['Runtime'] =  df['Cycles'] / (3.4*10**9)

data_perf = df.pivot_table('Performance', 'Block size X', 'Block size Y', fill_value=float("nan"))
data_runtime = df.pivot_table('Runtime', 'Block size X', 'Block size Y', fill_value=float("nan"))

print(data_runtime)

# Create the heatmap
x_size = 15/19 * data_perf.shape[1]
y_size = 10/17 * data_perf.shape[0]
plt.figure(figsize=(x_size, y_size)) # big plot
# plt.figure(figsize=(8, 5)) # smaller plot

if plot_type.startswith("perf"):
    sns.heatmap(data_perf, annot=True, fmt=".2f", cmap="YlOrRd")
else:
    # reverse colors cuz lower (runtime) is better
    sns.heatmap(data_runtime, annot=True, fmt=".2f", cmap="YlOrRd_r")

if plot_type.startswith("perf"):
    plt.title(f"Ya moma's performance (flop/cycle) for matrix_dim={matrix_dimension} & timestamp={timestamp}")
else:
    plt.title(f"Ya moma's runtime (sec) for matrix_dim={matrix_dimension} & timestamp={timestamp}")

base_dir = measurements_file.parent
if plot_type.startswith("perf"):
    heatmap_path = base_dir / 'heatmap_rect_perf.png'
else:
    heatmap_path = base_dir / 'heatmap_rect_runtime.png'

plt.savefig(heatmap_path, dpi=300)
print(f"Plot saved as '{heatmap_path}'")
