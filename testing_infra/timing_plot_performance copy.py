import matplotlib.pyplot as plt
import numpy as np

# Sample data based on the plot in the image
input_sizes = [16, 64, 256, 1_000, 4_000, 16_000, 64_000, 256_000, 1_000_000]
numerical_recipes = [0.5, 1, 1.2, 1.5, 1.7, 1.8, 1.9, 2, 2.1]
best_scalar_code = [1, 2, 3, 4, 4.5, 5, 5.2, 5.5, 5.7]
best_vector_code = [2, 4, 6, 8, 10, 12, 13, 14, 15]
best_vector_parallel_code = [3, 6, 12, 24, 30, 35, 32, 28, 25]

plt.figure(figsize=(10, 6))

# Plotting the lines and storing the line objects
line1, = plt.plot(input_sizes, numerical_recipes, marker='o', label='Numerical recipes')
line2, = plt.plot(input_sizes, best_scalar_code, marker='o', label='Best scalar code')
line3, = plt.plot(input_sizes, best_vector_code, marker='o', label='Best vector code')
line4, = plt.plot(input_sizes, best_vector_parallel_code, marker='o', label='Best code')

plt.xscale('log')

# Left-aligned title
plt.title(r'$\bf{DFT\ (single\ precision)\ on\ Intel\ Core\ i7\ (4\ cores)}$' + '\nPerformance [Flops/cycle] vs. input size', loc='left')


# Adding labels directly to the lines within the graph with matching colors
offset = 1.5  # Offset to move the text slightly above the lines
for label, y_values, line in zip(['Numerical recipes', 'Best scalar code', 'Best vector code', 'Best vector and\n parallel code'],
                                 [numerical_recipes, best_scalar_code, best_vector_code, best_vector_parallel_code],
                                 [line1, line2, line3, line4]):
    plt.text(input_sizes[6], y_values[6] + offset, label, fontsize=10, fontweight='bold', verticalalignment='bottom', color=line.get_color())

# Display grid
plt.grid(True, which='major', axis='y')


# Save the plot to a file
# plt.savefig("plot_perf.png", dpi=300)
# Save the plot to a file with specific dimensions
plt.savefig("plot_perf.pdf", format='pdf', bbox_inches='tight', pad_inches=0.1)


# Show plot
plt.show()
