import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Frequency of the processor
FREQUENCY = 2e9 

def plot_csv_files(csv_files, labels):
    # Set Seaborn style and context with light gray background
    sns.set_theme(style="dark")

    # Create a color palette
    palette = sns.color_palette("rocket", len(csv_files))

    for csv_file, label, color in zip(csv_files, labels, palette):
        # Read CSV file into a pandas DataFrame
        df = pd.read_csv(csv_file)

        df['Y'] /= FREQUENCY

        # Plot the data from the DataFrame with custom label
        sns.lineplot(data=df, x="X", y="Y", label=label, color=color)

    # Add labels, legend, and grid
    plt.xlabel("Sizes")
    plt.ylabel("Runtime [s]")
    plt.legend()
    plt.grid(axis='y', alpha=0.7) 

    # Show the plot
    plt.show()

# Example usage
csv_files = ["baseline.csv", "preallocated.csv"]  # List of CSV files
custom_labels = ["Baseline", "Preallocated"]  # List of custom labels
plot_csv_files(csv_files, custom_labels)
