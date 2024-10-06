import csv
import pandas as pd
import matplotlib.pyplot as plt
import sys
import os
import numpy as np

# Function to extract data from text file and append to a DataFrame
def extract_to_df(txt_file):
    data = {'Test': [], 'Size': [], 'Value': []}
    current_test = None
    with open(txt_file, 'r') as infile:
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
                    data['Test'].append(current_test)
                    data['Size'].append(float(size))
                    data['Value'].append(float(value))
    return pd.DataFrame(data)

# Ensure file paths are provided as command-line arguments
if len(sys.argv) < 2:
    print("Usage: python plot.py <path_to_txt_file1> <path_to_txt_file2> ...")
    sys.exit(1)

# Create necessary directories for iteration CSVs and aggregated CSVs
iterations_dir = 'csv_per_iterations'
aggregated_dir = 'csv_aggregated'
plots_dir = 'plots'

if not os.path.exists(iterations_dir):
    os.makedirs(iterations_dir)

if not os.path.exists(aggregated_dir):
    os.makedirs(aggregated_dir)

if not os.path.exists(plots_dir):
    os.makedirs(plots_dir)

# Dictionary to store dataframes for each case
case_data = {}

# Group files by case and read in data
for txt_file in sys.argv[1:]:
    # Extract the base name without the iteration part (i.e., everything before the last underscore and numbers)
    base_name = '_'.join(os.path.basename(txt_file).split('_')[:-1])  # Full name without iteration number
    
    # Extract data from text file to DataFrame
    df = extract_to_df(txt_file)
    
    # Save individual iteration CSV files
    iteration_csv_file = f"{iterations_dir}/{os.path.basename(txt_file)}.csv"
    df.to_csv(iteration_csv_file, index=False)
    
    if base_name not in case_data:
        case_data[base_name] = []
    
    # Append dataframe for the current iteration
    case_data[base_name].append(df)

# Initialize lists to store all latency and bandwidth data for plotting
latency_data = []
bandwidth_data = []
latency_std_values = []
bandwidth_std_values = []

# Process each case
for case, df_list in case_data.items():
    # Concatenate all iterations
    combined_df = pd.concat(df_list)
    
    # Calculate mean and std deviation for each size and test type
    grouped = combined_df.groupby(['Test', 'Size']).agg(
        Value_mean=('Value', 'mean'),
        Value_std=('Value', 'std')
    ).reset_index()
    
    # Round the mean and std deviation to 2 decimal places
    grouped['Value_mean'] = grouped['Value_mean'].round(2)
    grouped['Value_std'] = grouped['Value_std'].round(2)

    # Sort the DataFrame so that 'Latency' appears before 'Bandwidth'
    grouped['Test'] = pd.Categorical(grouped['Test'], categories=['Latency', 'Bandwidth'], ordered=True)
    grouped = grouped.sort_values(['Test', 'Size'])
    
    # Save aggregated CSV file for mean and std, using the full base name
    aggregated_csv_file = f"{aggregated_dir}/{case}_aggregated.csv"
    grouped.to_csv(aggregated_csv_file, index=False)
    
    # Separate data by test type for combined plotting
    latency_df = grouped[grouped['Test'] == 'Latency']
    bandwidth_df = grouped[grouped['Test'] == 'Bandwidth']
    
    if not latency_df.empty:
        latency_data.append((latency_df, case))
        latency_std_values.extend(latency_df['Value_std'].values)
    
    if not bandwidth_df.empty:
        bandwidth_data.append((bandwidth_df, case))
        bandwidth_std_values.extend(bandwidth_df['Value_std'].values)

# Color map for standard deviation visualization
cmap = plt.cm.viridis

# Different markers for each case
markers = ['o', 'x', 's']  # Circle, cross, square

# Get min and max std values for each test type
min_latency_std = min(latency_std_values)
max_latency_std = max(latency_std_values)
min_bandwidth_std = min(bandwidth_std_values)
max_bandwidth_std = max(bandwidth_std_values)

# Plot combined Latency data with raw std deviation color coding (specific to Latency)
plt.figure(figsize=(10, 6))
for i, (latency_df, label) in enumerate(latency_data):
    # Use raw std values for color scale specific to Latency test
    colors = cmap((latency_df['Value_std'] - min_latency_std) / (max_latency_std - min_latency_std))

    # Count the number of iterations for this case
    iteration_count = len(case_data[label])

    # Plot lines connecting the points
    plt.plot(latency_df['Size'], latency_df['Value_mean'], linestyle='-', color='gray', alpha=0.5)
    
    # Plot color-coded points based on std
    plt.scatter(latency_df['Size'], latency_df['Value_mean'], c=latency_df['Value_std'], cmap=cmap,
                vmin=min_latency_std, vmax=max_latency_std, marker=markers[i % len(markers)], 
                label=f'{label}', s=100, edgecolor='black')

# Update title to include iteration count
plt.xscale('log')
plt.yscale('log')
plt.title(f'Latency Test - Standard Deviation Visualization (Iterations: {iteration_count})')
plt.xlabel('Message Size (Bytes)')
plt.ylabel('Latency (us)')
plt.grid(True, which="both", ls="--")
plt.legend()
latency_plot_filename = f"{plots_dir}/latency_std_visualization.png"
plt.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=min_latency_std, vmax=max_latency_std)), 
             label='Standard Deviation (Latency)')
plt.savefig(latency_plot_filename)
plt.close()

# Plot combined Bandwidth data with raw std deviation color coding (specific to Bandwidth)
plt.figure(figsize=(10, 6))
for i, (bandwidth_df, label) in enumerate(bandwidth_data):
    # Use raw std values for color scale specific to Bandwidth test
    colors = cmap((bandwidth_df['Value_std'] - min_bandwidth_std) / (max_bandwidth_std - min_bandwidth_std))

    # Count the number of iterations for this case
    iteration_count = len(case_data[label])

    # Plot lines connecting the points
    plt.plot(bandwidth_df['Size'], bandwidth_df['Value_mean'], linestyle='-', color='gray', alpha=0.5)
    
    # Plot color-coded points based on std
    plt.scatter(bandwidth_df['Size'], bandwidth_df['Value_mean'], c=bandwidth_df['Value_std'], cmap=cmap,
                vmin=min_bandwidth_std, vmax=max_bandwidth_std, marker=markers[i % len(markers)], 
                label=f'{label}', s=100, edgecolor='black')

# Update title to include iteration count
plt.xscale('log')
plt.title(f'Bandwidth Test - Standard Deviation Visualization (Iterations: {iteration_count})')
plt.xlabel('Message Size (Bytes)')
plt.ylabel('Bandwidth (MB/s)')
plt.grid(True, which="both", ls="--")
plt.legend()
bandwidth_plot_filename = f"{plots_dir}/bandwidth_std_visualization.png"
plt.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=min_bandwidth_std, vmax=max_bandwidth_std)), 
             label='Standard Deviation (Bandwidth)')
plt.savefig(bandwidth_plot_filename)
plt.close()

print(f"CSV files for each iteration saved in: {iterations_dir}/")
print(f"Aggregated CSV files with mean and std saved in: {aggregated_dir}/")
print(f"Combined plots saved in: {plots_dir}/")