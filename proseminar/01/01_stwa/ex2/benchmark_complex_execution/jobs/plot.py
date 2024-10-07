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

# Ensure file paths are provided as command-line arguments
if len(sys.argv) < 2:
    print("Usage: python plot.py <path_to_txt_file1> <path_to_txt_file2> ...")
    sys.exit(1)

# Create 'plots' and 'csv' directories if they don't exist
plots_dir = 'plots'
csv_dir = 'csv'

if not os.path.exists(plots_dir):
    os.makedirs(plots_dir)

if not os.path.exists(csv_dir):
    os.makedirs(csv_dir)

# Initialize lists to store all latency and bandwidth data
latency_data = []
bandwidth_data = []

# Process each file provided as an argument
for txt_file in sys.argv[1:]:
    base_name = os.path.splitext(os.path.basename(txt_file))[0]
    csv_file = f"{csv_dir}/{base_name}.csv"
    
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

    # Append data to the lists with a label
    latency_data.append((latency_df, os.path.basename(txt_file)))
    bandwidth_data.append((bandwidth_df, os.path.basename(txt_file)))

# Plot Latency data from all files
plt.figure(figsize=(10, 6))
for latency_df, label in latency_data:
    plt.plot(latency_df['Size'], latency_df['Value'], marker='o', label=f'Latency - {label}')
plt.xscale('log')
plt.yscale('log')
plt.title('OSU MPI Latency Test')
plt.xlabel('Message Size (Bytes)')
plt.ylabel('Latency (us)')
plt.grid(True, which="both", ls="--")
plt.legend()
latency_plot_filename = f"{plots_dir}/latency_plot_all_cases.png"
plt.savefig(latency_plot_filename)
plt.close()

# Plot Bandwidth data from all files
plt.figure(figsize=(10, 6))
for bandwidth_df, label in bandwidth_data:
    plt.plot(bandwidth_df['Size'], bandwidth_df['Value'], marker='o', label=f'Bandwidth - {label}')
plt.xscale('log')
plt.title('OSU MPI Bandwidth Test')
plt.xlabel('Message Size (Bytes)')
plt.ylabel('Bandwidth (MB/s)')
plt.grid(True, which="both", ls="--")
plt.legend()
bandwidth_plot_filename = f"{plots_dir}/bandwidth_plot_all_cases.png"
plt.savefig(bandwidth_plot_filename)
plt.close()

print(f"Plots saved: \n{latency_plot_filename}\n{bandwidth_plot_filename}")
print(f"CSV files saved in: {csv_dir}/")
