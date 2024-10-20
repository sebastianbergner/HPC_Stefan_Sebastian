import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
import os

plt.rcParams['figure.dpi'] = 400

# Function to process the input CSV, compute metrics (mean, std, speedup, efficiency), and return processed data
def process_data(csv_file):
    # Read the input CSV file
    df = pd.read_csv(csv_file)

    # Aggregate by 'Impl/Ranks' and 'Problem Size' to calculate mean, standard deviation, and count repetitions
    aggregated_df = df.groupby(['Impl/Ranks', 'Problem Size'], sort=False).agg(
        Time_mean=('Time', 'mean'),
        Time_std=('Time', 'std'),
        Repetitions=('Time', 'size')  # Count the number of repetitions
    ).reset_index()

    # Preserve the original order of 'Impl/Ranks'
    impl_ranks_order = df['Impl/Ranks'].unique()
    aggregated_df['Impl/Ranks'] = pd.Categorical(aggregated_df['Impl/Ranks'], categories=impl_ranks_order, ordered=True)
    aggregated_df = aggregated_df.sort_values(by=['Impl/Ranks', 'Problem Size'])

    # Calculate rank from 'Impl/Ranks' (e.g., 'par_block/2' => 2)
    aggregated_df['Rank'] = aggregated_df['Impl/Ranks'].apply(lambda x: int(x.split('/')[-1]) if '/' in x else 1).astype(int)

    # Calculate speedup and efficiency using the seq/1 value for each problem size
    seq_df = aggregated_df[aggregated_df['Impl/Ranks'] == 'seq/1'][['Problem Size', 'Time_mean']].rename(columns={'Time_mean': 'Seq_Time'})
    merged_df = pd.merge(aggregated_df, seq_df, on='Problem Size')
    merged_df['Speedup'] = merged_df['Seq_Time'] / merged_df['Time_mean']
    merged_df['Efficiency'] = merged_df['Speedup'] / merged_df['Rank']

    # Round numeric columns to two decimal places for cleaner output
    merged_df = merged_df.round({'Time_mean': 2, 'Time_std': 2, 'Speedup': 2, 'Efficiency': 2})

    # Reorder the columns to place 'Repetitions' directly after 'Problem Size'
    merged_df = merged_df[['Impl/Ranks', 'Problem Size', 'Repetitions', 'Time_mean', 'Time_std', 'Speedup', 'Efficiency']]

    # Save the processed data to a new CSV file
    output_csv = f"{os.path.splitext(csv_file)[0]}_processed.csv"
    merged_df.to_csv(output_csv, index=False)
    print(f'Processed data saved to {output_csv}')

    # Return the processed data and maximum number of repetitions
    return output_csv, merged_df['Repetitions'].max(), merged_df

# Function to generate plots for Wall Time, Speedup, and Efficiency based on processed data
def create_plots(processed_csv, repetitions):
    # Read the processed data
    df = pd.read_csv(processed_csv)

    # Define color palette for the plots
    colors = ['#00429d', '#93c2ff', '#ff004d', '#ffa07a', '#f4e842', 
              '#3ebf72', '#00b7ff', '#ffb14e', '#ff7070', '#6a0dad', 
              '#ff77ff', '#008080', '#ff4500']

    # Create a Wall Time plot (log scale for time)
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.barplot(x='Impl/Ranks', y='Time_mean', hue='Problem Size', data=df, ax=ax, errorbar='sd', palette=colors)
    ax.set_yscale('log')
    ax.set_title(f'Wall Time (Mean) - {repetitions} Repetitions')
    ax.set_xlabel('Impl/Ranks')
    ax.set_ylabel('Time (log scale)')
    ax.tick_params(axis='x', rotation=45)
    plt.tight_layout()
    plt.savefig(f"{processed_csv.replace('_processed.csv', '_wall_time.png')}")
    plt.close()

    # Create a Speedup plot
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.barplot(x='Impl/Ranks', y='Speedup', hue='Problem Size', data=df, ax=ax, palette=colors)
    ax.set_title(f'Speedup - {repetitions} Repetitions')
    ax.set_xlabel('Impl/Ranks')
    ax.set_ylabel('Speedup')
    ax.tick_params(axis='x', rotation=45)
    plt.tight_layout()
    plt.savefig(f"{processed_csv.replace('_processed.csv', '_speedup.png')}")
    plt.close()

    # Create an Efficiency plot
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.barplot(x='Impl/Ranks', y='Efficiency', hue='Problem Size', data=df, ax=ax, palette=colors)
    ax.set_title(f'Efficiency - {repetitions} Repetitions')
    ax.set_xlabel('Impl/Ranks')
    ax.set_ylabel('Efficiency')
    ax.tick_params(axis='x', rotation=45)
    plt.tight_layout()
    plt.savefig(f"{processed_csv.replace('_processed.csv', '_efficiency.png')}")
    plt.close()

# Function to generate a LaTeX table based on processed data in the required format
def create_latex_table(df, latex_output_file):
    # Get all unique problem sizes
    unique_problem_sizes = sorted(df['Problem Size'].unique())

    with open(latex_output_file, 'w') as f:
        # Begin LaTeX table
        f.write("\\begin{tabular}{|" + "l" + "l" * (4 * len(unique_problem_sizes)) + "|}\n")
        f.write("\\hline\n")
        f.write("\\multicolumn{" + str(4 * len(unique_problem_sizes) + 1) + "}{|c|}{\\textbf{Results of 1D Heat Stencil Execution}} \\\\ \\hline\n")
        
        # First header row with "Impl/Ranks" and "Problem Size"
        f.write("\\multicolumn{1}{|c|}{\\textbf{Impl/Ranks}} & " +
                "\\multicolumn{" + str(4 * len(unique_problem_sizes)) + "}{c|}{\\textbf{Problem Size}} \\\\ \\hline\n")
        
        # Second header row with Problem Sizes (in 4 columns per size)
        f.write("\\multicolumn{1}{|c|}{\\textbf{}} & " +
                " & ".join([f"\\multicolumn{{4}}{{c|}}{{\\textbf{{{int(psize)}}}}}" for psize in unique_problem_sizes]) +
                " \\\\ \\hline\n")
        
        # Third header row with μ, σ, S(p), E(p)
        f.write("\\multicolumn{1}{|l|}{} & " +
                " & ".join(["\\multicolumn{1}{c|}{$\\mu$ [s]} & \\multicolumn{1}{c|}{$\\sigma$ [s]} & \\multicolumn{1}{c|}{S(p)} & \\multicolumn{1}{c|}{E(p)}"
                            for _ in unique_problem_sizes]) +
                " \\\\ \\hline\n")

        # Iterate through each unique Impl/Ranks and write corresponding rows
        for impl_rank in df['Impl/Ranks'].unique():
            impl_rank_escaped = impl_rank.replace('_', '\\_')  # Escape underscores in LaTeX
            f.write(f"\\multicolumn{{1}}{{|l|}}{{{impl_rank_escaped}}} ")
            for psize in unique_problem_sizes:
                # Get the row data for the given Impl/Ranks and Problem Size
                row = df[(df['Impl/Ranks'] == impl_rank) & (df['Problem Size'] == psize)]
                if not row.empty:
                    mean = row['Time_mean'].values[0]
                    std = row['Time_std'].values[0]
                    speedup = row['Speedup'].values[0]
                    efficiency = row['Efficiency'].values[0]
                    f.write(f" & \\multicolumn{{1}}{{r|}}{{{mean:.2f}}} & \\multicolumn{{1}}{{r|}}{{{std:.2f}}} & \\multicolumn{{1}}{{r|}}{{{speedup:.2f}}} & \\multicolumn{{1}}{{r|}}{{{efficiency:.2f}}} ")
                else:
                    f.write(" & \\multicolumn{1}{l|}{-} & \\multicolumn{1}{l|}{-} & \\multicolumn{1}{l|}{-} & \\multicolumn{1}{l|}{-} ")  # Use dashes if no data available
            f.write(" \\\\ \\hline\n")  # End row with horizontal line

        # End LaTeX table
        f.write("\\end{tabular}\n")

    print(f'LaTeX table saved to {latex_output_file}')


# Main function to process CSV, generate plots, and create a LaTeX table
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process data, generate plots, and create a LaTeX table.')
    parser.add_argument('csv', nargs='*', type=str, help='List of CSV files to process and generate output from')
    args = parser.parse_args()

    csv_files = args.csv
    for csv_file in csv_files:
        # Step 1: Process the CSV data
        processed_csv, repetitions, processed_data = process_data(csv_file)

        # Step 2: Generate plots for Wall Time, Speedup, and Efficiency
        create_plots(processed_csv, repetitions)

        # Step 3: Create a LaTeX table based on the processed data
        latex_output = processed_csv.replace('_processed.csv', '_table.tex')
        create_latex_table(processed_data, latex_output)
