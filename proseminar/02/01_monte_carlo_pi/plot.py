import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import argparse

plt.rcParams['figure.dpi'] = 400


def sns_plotter(x, y, ax, title, xlabel, ylabel, df, hue, errorbar='sd'):
    plot = sns.barplot(x=x, y=y, hue=hue, data=df, ax=ax, errorbar=errorbar)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.tick_params(axis='x', rotation=0)
    return plot

def plot_all(filename: str):
    # read data
    # Read the CSV file into a DataFrame
    df = pd.read_csv(filename)

    df['Problem size'] = pd.to_numeric(df['Problem size'], errors='coerce')
    df['Time'] = pd.to_numeric(df['Time'], errors='coerce')

    # Bar Plots for Mean values
    fig, axs = plt.subplots(1, 1, figsize=(12, 10))

    # Wall Time Mean
    sns_plotter(x='Type', y='Time', hue='Problem size', df=df, ax=axs, errorbar='sd', title='Measurement', xlabel='Type', ylabel='Time (log)')

    plt.yscale('log')

    plt.grid(True, which='both', axis='y', linestyle='-', linewidth=0.3)
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