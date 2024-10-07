import csv
import pandas as pd
import matplotlib.pyplot as plt
import sys
import os

# Function to extract data from text file and write to CSV
def extract_to_csv(txt_file, csv_file):
    with open(txt_file, 'r') as infile, open(csv_file, 'w', newline='') as outfile:
        writer = csv.writer(outfile)
        writer.writerow(['Test', 'Size', 'Value'])  # Writing header
        current_test = None
        
        for line in infile:
            line = line.strip()
            if "Latency" in line:
                current_test = 'Latency'
            elif "Bandwidth" in line:
                current_test = 'Bandwidth'
            elif line and not line.startswith('#') and not line.startswith("mpiexec"):
                parts = line.split()
                if len(parts) == 2 and current_test:
                    size, value = parts
                    writer.writerow([current_test, size, value])

# Ensure a file path is provided as a command-line argument
if len(sys.argv) != 2:
    print("Usage: python plot.py <path_to_txt_file>")
    sys.exit(1)

# Get the file path from the command line argument
txt_file = sys.argv[1]

# Create a CSV file path with the same base name as the txt file
base_name = os.path.splitext(txt_file)[0]
csv_file = f"{base_name}.csv"

# Extract data from text file to CSV
extract_to_csv(txt_file, csv_file)

# Read the CSV file into a DataFrame
df = pd.read_csv(csv_file)

# Convert Size and Value columns to numeric and clean the data
df['Size'] = pd.to_numeric(df['Size'], errors='coerce')
df['Value'] = pd.to_numeric(df['Value'], errors='coerce')
df = df.dropna()

# Separate the data by test type
latency_df = df[df['Test'] == 'Latency']
bandwidth_df = df[df['Test'] == 'Bandwidth']

# Create 'plots' directory if it doesn't exist
plots_dir = 'plots'
if not os.path.exists(plots_dir):
    os.makedirs(plots_dir)

# Plot for Latency Test
plt.figure(figsize=(10, 6))
plt.plot(latency_df['Size'], latency_df['Value'], marker='o', label='Latency (us)')
plt.xscale('log')
plt.yscale('log')
plt.title('OSU MPI Latency Test')
plt.xlabel('Message Size (Bytes)')
plt.ylabel('Latency (us)')
plt.grid(True, which="both", ls="--")
plt.legend()
latency_plot_filename = f"{plots_dir}/latency_plot.png"
plt.savefig(latency_plot_filename)  # Save the plot as a file
plt.close()  # Close the plot to avoid showing it

# Plot for Bandwidth Test
plt.figure(figsize=(10, 6))
plt.plot(bandwidth_df['Size'], bandwidth_df['Value'], marker='o', label='Bandwidth (MB/s)', color='orange')
plt.xscale('log')
plt.title('OSU MPI Bandwidth Test')
plt.xlabel('Message Size (Bytes)')
plt.ylabel('Bandwidth (MB/s)')
plt.grid(True, which="both", ls="--")
plt.legend()
bandwidth_plot_filename = f"{plots_dir}/bandwidth_plot.png"
plt.savefig(bandwidth_plot_filename)  # Save the plot as a file
plt.close()  # Close the plot to avoid showing it

print(f"Plots saved: \n{latency_plot_filename}\n{bandwidth_plot_filename}")
