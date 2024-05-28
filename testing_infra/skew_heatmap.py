
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

df = pd.read_csv('performance_metrics.csv')

df['GFLOP Rate'] = df['FLOPs'] / df['Cycles']

# Dimensions for the x and y axes
x_labels = [4, 8, 16, 32, 64, 128, 256]
y_labels = [4, 8, 16, 32, 64, 128, 256]

# Create the heatmap
plt.figure(figsize=(8, 6))
sns.heatmap(data, annot=True, fmt=".1f", cmap="YlOrRd", xticklabels=x_labels, yticklabels=y_labels)

# Add labels and title
plt.xlabel("X-Dimension of Cache Block")
plt.ylabel("Y-Dimension of Cache Block")
plt.title("Iteration #1: GFlop Rate")

# Show the heatmap
plt.show()
plt.savefig("heatmap", dpi=300)
