import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Frequency of the processor
FREQUENCY = 2e9 

def plot_csv_files(csv_files, labels):

    sns.set_theme(style="dark")

    palette = sns.color_palette("rocket", len(csv_files))

    for csv_file, label, color in zip(csv_files, labels, palette):
        df = pd.read_csv(csv_file)
        df['Y'] /= FREQUENCY
        sns.lineplot(data=df, x="X", y="Y", label=label, color=color)

    plt.xlabel("Sizes")
    plt.ylabel("Runtime [s]")
    plt.legend()
    plt.grid(axis='y', alpha=0.7) 


    plt.show()


csv_files = ["baseline.csv", "preallocated.csv"]  # List of CSV files
custom_labels = ["Baseline", "Preallocated"]  # List of custom labels
plot_csv_files(csv_files, custom_labels)
