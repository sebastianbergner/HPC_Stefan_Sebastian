import pandas as pd
import argparse
import os

# Function to aggregate, calculate std, round values, and preserve order of 'Impl/Ranks'
def aggregate_and_calculate_std(input_csv):
    # Read the CSV file
    df = pd.read_csv(input_csv)

    # Group by 'Impl/Ranks' and 'Problem Size', then calculate mean and std
    aggregated_df = df.groupby(['Impl/Ranks', 'Problem Size'], sort=False).agg(
        Time_mean=('Time', 'mean'),
        Time_std=('Time', 'std')
    ).reset_index()

    # Maintain the order of 'Impl/Ranks' as it appears in the input file
    impl_ranks_order = df['Impl/Ranks'].unique()
    aggregated_df['Impl/Ranks'] = pd.Categorical(aggregated_df['Impl/Ranks'], categories=impl_ranks_order, ordered=True)
    aggregated_df = aggregated_df.sort_values(by=['Impl/Ranks', 'Problem Size'])

    # Calculate rank from 'Impl/Ranks'
    aggregated_df['Rank'] = aggregated_df['Impl/Ranks'].apply(lambda x: int(x.split('/')[-1]) if '/' in x else 1).astype(int)
    
    # Create a temporary dataframe for seq/1 times for each problem size
    seq_df = aggregated_df[aggregated_df['Impl/Ranks'] == 'seq/1'][['Problem Size', 'Time_mean']].rename(columns={'Time_mean': 'Seq_Time'})
    
    # Merge seq/1 times with the aggregated dataframe
    merged_df = pd.merge(aggregated_df, seq_df, on='Problem Size')
    
    # Calculate speedup and efficiency
    merged_df['Speedup'] = merged_df['Seq_Time'] / merged_df['Time_mean']
    merged_df['Efficiency'] = merged_df['Speedup'] / merged_df['Rank']
    
    # Round all numeric columns to 2 decimal places
    merged_df = merged_df.round({'Time_mean': 2, 'Time_std': 2, 'Speedup': 2, 'Efficiency': 2})
    
    # Drop Seq_Time column (optional)
    merged_df = merged_df.drop(columns=['Seq_Time'])

    # Drop Ranks if you don't want them in the final output
    merged_df = merged_df.drop(columns=['Rank'])

    # Generate the output filename by appending '_aggregated.csv'
    base_name = os.path.splitext(input_csv)[0]
    output_csv = f"{base_name}_aggregated.csv"

    # Save the results to the new CSV file
    merged_df.to_csv(output_csv, index=False)

    print(f'Aggregated results saved to {output_csv}')

# Main function to handle command line arguments
def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Aggregate results, calculate std, speedup, efficiency, and preserve the order of Impl/Ranks')
    
    # Define arguments
    parser.add_argument('input_csv', help='Path to the input CSV file')

    # Parse arguments
    args = parser.parse_args()

    # Call the function with parsed arguments
    aggregate_and_calculate_std(args.input_csv)

# Entry point
if __name__ == '__main__':
    main()
