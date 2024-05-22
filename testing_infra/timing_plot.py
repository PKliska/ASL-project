import sys

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

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
        df['n_cycles'] = df['n_cycles'] / FREQUENCY
        sns.lineplot(data=df, x="matrix_dimension", y="n_cycles", label=label)
        # sns.lineplot(data=df, x="matrix_dimension", y="n_cycles", label=label, color=color)

    plt.xlabel("Dimension")
    plt.ylabel("Runtime")
    plt.legend()
    plt.grid(axis='y', alpha=0.7)

    # add show ticks for matrix_dimension for which we run the algo
    plt.xticks(df["matrix_dimension"])

    # add vertical lines at these places so that it's easier to see
    for x in df["matrix_dimension"]:
        plt.axvline(x=x, color='w', linestyle='-', linewidth=1)

    # plt.figure(figsize=(10, 6))
    plt.savefig("plot", dpi=300)
    plt.show()

csv_files = (sys.argv[1]).split(",")  # List of CSV files, comma separeted
custom_labels = [filename.removesuffix(".csv") for filename in csv_files]  # List of custom labels
plot_csv_files(csv_files, custom_labels)

