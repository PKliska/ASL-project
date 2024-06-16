import sys

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import os

# Frequency of the processor
FREQUENCY = 3.4e9 # 3.4 GHz

def plot_csv_files(csv_files, labels, label_indices):
    # sns.set_theme(style="d¨ark")

    palette = sns.color_palette("deep", len(csv_files))
    lines = []

    for csv_file, label, index, color in zip(csv_files, labels, label_indices, palette):
        # flops (91 * (nx-2)*(ny-2) + 22 * (nx-2)(ny-2)*pit + 4) * steps
        # for us nx=ny=41 and pit = 50, so
        # num_flops = 1811515 * steps
        df = pd.read_csv(csv_file)
        #performance = 1811515 * df['matrix_dimension'] / df['n_cycles']
        # flops = 76_000_000_000

        df['n_performance'] = df['n_flops'] / df['n_cycles']
        # to show runtime (sec) it's the below:
        # df['n_cycles'] = df['n_cycles'] / FREQUENCY

        # sns.lineplot(data=df, x="matrix_dimension", y="n_performance", label=label)
        # sns.lineplot(data=df, x="matrix_dimension", y="n_cycles", label=label, color=color)
         # Plotting the line and storing the line object
        line, = plt.plot(df['matrix_dimension'], df['n_performance'], marker='o',  markersize=3, label=label, color=color)
        lines.append((label, df['n_performance'], line, index))
    
       # Adding labels directly to the lines within the graph with matching colors
    offset = 0.05 # Offset to move the text slightly above the lines
    for label, y_values, line, index in lines:
        plt.text(df['matrix_dimension'].iloc[index], y_values.iloc[index] + offset, label, fontsize=10, fontweight='bold', verticalalignment='bottom', color=line.get_color())

    plt.title('Performance [Flops/cycle] vs. input size', loc='left')
    plt.grid(axis='y', alpha=0.7)

    # l1 = np.sqrt(32*2**10/8)
    # l2 = np.sqrt(256*2**10/8)
    # l3 = np.sqrt(8*2**20/8)
    # plt.vlines(l1, ymin=0.5, ymax=3.5, color='red')
    # plt.annotate('L1', (l1, 1.0), color='red', xytext=(5,0), textcoords='offset pixels')
    # plt.vlines(l2, ymin=0.5, ymax=3.5, color='red')
    # plt.annotate('L2', (l2, 1.0), color='red', xytext=(5,0), textcoords='offset pixels')
    # plt.vlines(l3, ymin=0.5, ymax=3.5, color='red')
    # plt.annotate('L3', (l3, 1.0), color='red', xytext=(5,0), textcoords='offset pixels')


    # add show ticks for matrix_dimension for which we run the algo
    # selected_indices = [0, 4, 6, 8, 11]  # 0-based indices
    # selected_matrix_dimensions = df['matrix_dimension'].iloc[selected_indices]
    # plt.xticks(ticks=selected_matrix_dimensions, labels=selected_matrix_dimensions)
    xticks_indices = list(range(0, len(df['matrix_dimension']), 2))
    if 7 not in xticks_indices:
        xticks_indices.append(2)
    xticks_indices.sort()

    plt.xticks(df['matrix_dimension'].iloc[xticks_indices])

    # plt.xticks(df['matrix_dimension'])

    # add vertical lines at these places so that it's easier to see
    # for x in df["matrix_dimension"]:
    #     plt.axvline(x=x, color='w', linestyle='-', linewidth=1)

    # plt.figure(figsize=(10, 6))
    plt.savefig("plot_perf_trapeze_icc_report", dpi=300)
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python3 timing_plot_performance_report.py <csv_files> <indices>")
        sys.exit(1)

    csv_files = sys.argv[1].split(",")  # List of CSV files, comma separated
    custom_labels = [filename.split("/")[-1].removesuffix(".csv") for filename in csv_files]  # Extract custom labels
    label_indices = list(map(int, sys.argv[2].split(",")))  # List of indices, comma separated

    if len(csv_files) != len(label_indices):
        print("Error: The number of indices must match the number of CSV files.")
        sys.exit(1)

    plot_csv_files(csv_files, custom_labels, label_indices)

