import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Frequency of the processor
FREQUENCY = 2.5e9 # 2.5 GHz

def plot_csv_files(csv_files, labels):
    sns.set_theme(style="dark")

    palette = sns.color_palette("rocket", len(csv_files))

    for csv_file, label, color in zip(csv_files, labels, palette):
        df = pd.read_csv(csv_file)
        df['n_cycles'] /= FREQUENCY
        sns.lineplot(data=df, x="n_simulation_iterations", y="n_cycles", label=label, color=color)

    plt.xlabel("n_simulation_iterations")
    plt.ylabel("Runtime [s]")
    plt.legend()
    plt.grid(axis='y', alpha=0.7)

    # add show ticks for n_simulation_iterations for which we run the algo
    plt.xticks(df["n_simulation_iterations"])

    # add vertical lines at these places so that it's easier to see
    for x in df["n_simulation_iterations"]:
        plt.axvline(x=x, color='w', linestyle='-', linewidth=1)

    plt.savefig("plot")
    # plt.show()

csv_files = ["baseline.csv", "preallocated.csv"]  # List of CSV files
custom_labels = ["Baseline", "Preallocated"]  # List of custom labels
plot_csv_files(csv_files, custom_labels)
