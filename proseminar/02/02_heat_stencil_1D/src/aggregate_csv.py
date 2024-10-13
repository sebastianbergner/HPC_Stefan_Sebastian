import pandas as pd
import argparse
import os

# Function to aggregate, calculate std, and sort alphabetically by 'Impl/Ranks'
def aggregate_and_calculate_std(input_csv):
    # Read the CSV file
    df = pd.read_csv(input_csv)

    # Group by 'Impl/Ranks' and 'Problem Size', then calculate mean and std
    aggregated_df = df.groupby(['Impl/Ranks', 'Problem Size']).agg(
        Time_mean=('Time', 'mean'),
        Time_std=('Time', 'std')
    ).reset_index()

    # Round the mean and std to two decimal places
    aggregated_df['Time_mean'] = aggregated_df['Time_mean'].round(2)
    aggregated_df['Time_std'] = aggregated_df['Time_std'].round(2)

    # Sort alphabetically by 'Impl/Ranks'
    aggregated_df = aggregated_df.sort_values(by='Impl/Ranks')

    # Generate the output filename by appending '_aggregated.csv'
    base_name = os.path.splitext(input_csv)[0]
    output_csv = f"{base_name}_aggregated.csv"

    # Save the results to the new CSV file
    aggregated_df.to_csv(output_csv, index=False)

    print(f'Aggregated results saved to {output_csv}')

# Main function to handle command line arguments
def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Aggregate results, calculate std, and sort alphabetically by Impl/Ranks')
    
    # Define arguments
    parser.add_argument('input_csv', help='Path to the input CSV file')

    # Parse arguments
    args = parser.parse_args()

    # Call the function with parsed arguments
    aggregate_and_calculate_std(args.input_csv)

# Entry point
if __name__ == '__main__':
    main()

