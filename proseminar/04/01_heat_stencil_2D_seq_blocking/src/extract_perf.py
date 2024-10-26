import re
import argparse


def main(input_file, output_file):
    # Regular expressions to extract the relevant information
    program_name_re = re.compile(r"for '(.+)'")
    metrics_re = {
        "L1_dcache_load_misses": re.compile(r"([\d,.]+)\s+L1-dcache-load-misses:u"),
        "L1_dcache_loads": re.compile(r"([\d,.]+)\s+L1-dcache-loads:u"),
        "LLC_load_misses": re.compile(r"([\d,.]+)\s+LLC-load-misses:u"),
        "LLC_loads": re.compile(r"([\d,.]+)\s+LLC-loads:u"),
        "time_elapsed": re.compile(r"([\d,.]+) seconds time elapsed"),
        "user_time": re.compile(r"([\d,.]+) seconds user"),
        "sys_time": re.compile(r"([\d,.]+) seconds sys")
    }

    # List to store the extracted data
    data = []

    # Parse the file and extract the metrics
    with open(input_file, 'r') as file:
        current_data = {}
        for line in file:
            # Detect program name
            program_match = program_name_re.search(line)
            if program_match:
                if current_data:  # if data already exists, append to list
                    data.append(current_data)
                    current_data = {}
                current_data['Program'] = program_match.group(1)

            # Extract each metric using regex
            for metric, regex in metrics_re.items():
                match = regex.search(line)
                if match:
                    current_data[metric] = match.group(1).replace('.', ':').replace(',', '.').replace(':', ',')

        if current_data:
            data.append(current_data)  # Append the last parsed program data

    # Generate LaTeX table
    with open(f"{output_file}.tex", 'w') as f:
        f.write("\\begin{table}[h!]\n")
        f.write("\\footnotesize\n")
        f.write("\\centering\n")
        f.write("\\begin{tabular}{|l|r|r|r|r|r|r|r|}\n")
        f.write("\\hline\n")
        f.write("Program & L1d load misses & L1d loads & LLC load misses & LLC loads & Time (s) & User time (s) & Sys time (s) \\\\ \\hline\n")

        for entry in data:
            f.write(f"{entry.get('Program', 'N/A')} & ")
            f.write(f"{entry.get('L1_dcache_load_misses', 'N/A')} & ")
            f.write(f"{entry.get('L1_dcache_loads', 'N/A')} & ")
            f.write(f"{entry.get('LLC_load_misses', 'N/A')} & ")
            f.write(f"{entry.get('LLC_loads', 'N/A')} & ")
            f.write(f"{entry.get('time_elapsed', 'N/A')} & ")
            f.write(f"{entry.get('user_time', 'N/A')} & ")
            f.write(f"{entry.get('sys_time', 'N/A')} \\\\ \\hline\n")

        f.write("\\end{tabular}\n")
        f.write("\\caption{Performance Metrics for Various Programs}\n")
        f.write("\\label{table:performance_metrics}\n")
        f.write("\\end{table}\n")

if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Process logs and create LaTeX tables.')
    parser.add_argument('logs', nargs='*', type=str, help='List of perf files to process and generate output from')
    args = parser.parse_args()

    files = args.logs
    for file in files:
        main(file, file.split('.')[0])