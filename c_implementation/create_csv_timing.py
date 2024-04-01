import subprocess
import matplotlib.pyplot as plt
import csv

# Define the compiler flags
compiler_flags = ["-I", "preallocated", "-t"]  # Example flags, modify as needed

# Run the compiled C program with compiler flags and capture its output
command = ["/Users/peterhorcic/team22/c_implementation/build/bin/cavity_flow"] + compiler_flags
process = subprocess.Popen(command, stdout=subprocess.PIPE)
output, _ = process.communicate()

# Decode the output and extract relevant data
output_lines = output.decode().split("\n")
# Assuming each line contains two numbers separated by space
data = [tuple(map(float, line.split())) for line in output_lines if line]

# Save data to a CSV file
with open("output.csv", "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["X", "Y"])  # Write header row
    writer.writerows(data)       # Write data rows

print("Data saved to output.csv")
