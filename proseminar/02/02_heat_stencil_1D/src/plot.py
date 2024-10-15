import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import argparse

plt.rcParams['figure.dpi'] = 400

def sns_plotter(x, y, ax, title, xlabel, ylabel, df, hue, errorbar='sd', palette=None):
    plot = sns.barplot(x=x, y=y, hue=hue, data=df, ax=ax, errorbar=errorbar, palette=palette)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.tick_params(axis='x', rotation=45)
    return plot

def plot_all(filename: str):
    # read data
    # Read the CSV file into a DataFrame
    df = pd.read_csv(filename)

    df['Problem Size'] = pd.to_numeric(df['Problem Size'], errors='coerce')
    df['Time'] = pd.to_numeric(df['Time'], errors='coerce')

    # Bar Plots for Mean values
    fig, axs = plt.subplots(1, 1, figsize=(12, 10))

    # Define a color palette for the bars
    colors = ['#00429d', '#93c2ff', '#ff004d', '#ffa07a', '#f4e842', 
          '#3ebf72', '#00b7ff', '#ffb14e', '#ff7070', '#6a0dad', 
          '#ff77ff', '#008080', '#ff4500']



    # Wall Time Mean
    sns_plotter(x='Impl/Ranks', y='Time', hue='Problem Size', df=df, ax=axs, errorbar='sd', 
                title='Measurement', xlabel='Implementation', ylabel='Time (log)', palette=colors)

    plt.yscale('log')

    plt.grid(True, which='both', axis='y', linestyle='-', linewidth=0.2)
    plt.gca().set_axisbelow(True)

    plt.tight_layout()
    plt.savefig(f'./{filename.replace(".csv", "")}_plot.png')
    plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot raw.csv files from exec script.')
    parser.add_argument('csv', nargs='*', type=str, help='name of all csv files')
    args = parser.parse_args()

    csv_files = args.csv
    for file in csv_files:
        plot_all(file)
