import subprocess
import matplotlib.pyplot as plt
import csv

# Define the compiler flags as desired
compiler_flags = ["-I", "preallocated", "-t"]  

# Modify this as required
command = ["/Users/peterhorcic/team22/c_implementation/build/bin/cavity_flow"] + compiler_flags
process = subprocess.Popen(command, stdout=subprocess.PIPE)
output, _ = process.communicate()

# Decode the output and extract relevant data
output_lines = output.decode().split("\n")
# Assuming each line contains two numbers separated by space
data = [tuple(map(float, line.split())) for line in output_lines if line]

with open("output.csv", "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["X", "Y"])  
    writer.writerows(data)       

print("Data saved to output.csv")
