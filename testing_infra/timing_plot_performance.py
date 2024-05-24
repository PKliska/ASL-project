import sys

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# Frequency of the processor
FREQUENCY = 3.4e9 # 3.4 GHz

def plot_csv_files(csv_files, labels):
    # sns.set_theme(style="dark")

    palette = sns.color_palette("rocket", len(csv_files))

    for csv_file, label, color in zip(csv_files, labels, palette):
        # flops (91 * (nx-2)*(ny-2) + 22 * (nx-2)(ny-2)*pit + 4) * steps
        # for us nx=ny=41 and pit = 50, so
        # num_flops = 1811515 * steps
        df = pd.read_csv(csv_file)
        #performance = 1811515 * df['matrix_dimension'] / df['n_cycles']
        # flops = 76_000_000_000

        df['n_performance'] = df['n_flops'] / df['n_cycles']
        # to show runtime (sec) it's the below:
        # df['n_cycles'] = df['n_cycles'] / FREQUENCY

        sns.lineplot(data=df, x="matrix_dimension", y="n_performance", label=label)
        # sns.lineplot(data=df, x="matrix_dimension", y="n_cycles", label=label, color=color)

    plt.xlabel("Dimension")
    plt.ylabel("Performance")
    plt.legend()
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
    plt.xticks(df["matrix_dimension"])

    # add vertical lines at these places so that it's easier to see
    for x in df["matrix_dimension"]:
        plt.axvline(x=x, color='w', linestyle='-', linewidth=1)

    # plt.figure(figsize=(10, 6))
    plt.savefig("plot_perf", dpi=300)
    plt.show()

csv_files = (sys.argv[1]).split(",")  # List of CSV files, comma separeted
custom_labels = [filename.removesuffix(".csv") for filename in csv_files]  # List of custom labels
plot_csv_files(csv_files, custom_labels)

