import pandas as pd
import matplotlib.pyplot as plt
import os
import sys

def plot_sequence(directory_path, sequence_number):
    """
    Plots the A_data and B_data of the sequence as bar plots with a common x-axis.
    Titles are positioned vertically on the right edge of the plots.

    Parameters:
    directory_path (str): The path to the directory containing the sequence files.
    sequence_number (int): The sequence number to plot.
    """
    sequence_file = os.path.join(directory_path, f"sequence_{sequence_number}.csv")

    if not os.path.exists(sequence_file):
        print(f"File not found: {sequence_file}")
        return

    # Load the data
    df = pd.read_csv(sequence_file)

    # Plotting
    fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(12, 6), sharex=True)

    # Plot A_data
    axes[0].bar(df.index, df['A_data'], color='blue')
    axes[0].set_ylabel('Value')

    # Plot B_data
    axes[1].bar(df.index, df['B_data'], color='red')
    axes[1].set_xlabel('Index')
    axes[1].set_ylabel('Value')

    # Set titles on the right edge of the plots
    axes[0].text(1.02, 0.5, 'A_data', va='center', ha='left', rotation='vertical', transform=axes[0].transAxes)
    axes[1].text(1.02, 0.5, 'B_data', va='center', ha='left', rotation='vertical', transform=axes[1].transAxes)

    # Adjust layout
    plt.subplots_adjust(hspace=0)
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python visualise_sequence.py <directory_path> <sequence_number>")
    else:
        directory_path = sys.argv[1]
        sequence_number = int(sys.argv[2])
        plot_sequence(directory_path, sequence_number)
