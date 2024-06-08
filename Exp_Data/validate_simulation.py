import numpy as np
import pandas as pd
import os
import sys
from multiprocessing import Pool
from contextlib import closing
import plotly.graph_objects as go
import matplotlib.pyplot as plt
from matplotlib import rcParams
import seaborn as sns

from validate_data import decode_with_viterbi, decode_with_k, decode_with_noanta, decode_with_k_patch, plot_sequences

# Ensure the parent directory is in the path for importing the Dataframe class
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(parent_dir)
from Epigenetic_Sequence_Predictor.ESP_Sim_Viterbi import Dataframe


def correct_and_conclude(df, correction_method, method_name, sim_path):
    """
    Corrects sequences using the specified method, generates conclusions,
    and saves the necessary files including plots.

    Parameters:
    - df: Dataframe instance containing sequences.
    - correction_method: Function to correct sequences.
    - method_name: Name of the correction method used.
    - sim_path: Path where the results and plots will be saved.
    """
    corrected_sequences = []

    # Correct each corrupt daughter sequence using the specified method
    for corrupt_daughter in df.corrupt_daughter_list:
        corrected_daughter = correction_method(df.alpha, df.beta, df.rho, df.mu, corrupt_daughter)
        corrected_sequences.append(corrected_daughter)

    # Update the Dataframe with the corrected sequences
    df.corrected_daughter_list = np.array(corrected_sequences)

    method_path = os.path.join(sim_path, method_name)

    # Save the conclusions and plot the sequences
    df.conclusions(method_path, level=0)
    # Note: Ensure plot_sequences is properly defined to accept sim_path or adjust as needed
    plot_sequences(method_path)

    # Return average BER for this method
    average_ber = np.mean(df.biterror)
    return method_name, average_ber


def plot_bit_error(csv_path):
    df = pd.read_csv(csv_path)
    plot_dir = os.path.join(os.path.dirname(csv_path), "plots")
    os.makedirs(plot_dir, exist_ok=True)
    methods = df['method'].unique()

    for method in methods:
        # Matplotlib Plot
        plt.figure(figsize=(7.5, 5.5), dpi=300)
        grouped = df[df['method'] == method].groupby(['alpha', 'beta'])
        for (alpha, beta), group in grouped:
            sorted_group = group.sort_values('mu')
            plt.plot(sorted_group['mu'].values, sorted_group['average_ber'].values, '-o',
                     label=f'Alpha={alpha}, Beta={beta}')
        plt.title(f'Method: {method} - Average BER vs $\mu$')
        plt.xlabel(r'$\mu$')
        plt.ylabel('Average BER')
        plt.legend(loc='best')
        plt.grid(True)
        for fmt in ['png', 'tiff', 'eps']:
            matplotlib_path = os.path.join(plot_dir, f'{method}_bit_error_vs_mu.{fmt}')
            plt.savefig(matplotlib_path, dpi=300, format=fmt)
        plt.close()
        # print(f'Matplotlib plot saved: {matplotlib_path}')

        # Plotly Plot
        fig = go.Figure()
        for (alpha, beta), group in grouped:
            sorted_group = group.sort_values('mu')
            fig.add_trace(go.Scatter(x=sorted_group['mu'], y=sorted_group['average_ber'], mode='lines+markers',
                                     name=f'Alpha={alpha}, Beta={beta}'))
        fig.update_layout(title=f'Method: {method} - Average BER vs $\mu$', xaxis_title='Mu',
                          yaxis_title='Average BER', legend_title="Alpha, Beta")
        plotly_path = os.path.join(plot_dir, f'{method}_bit_error_vs_mu.html')
        fig.write_html(plotly_path)
        # print(f'Plotly plot saved: {plotly_path}')


def plot_delta_error(csv_path):
    df = pd.read_csv(csv_path)
    plot_dir = os.path.join(os.path.dirname(csv_path), "plots")
    os.makedirs(plot_dir, exist_ok=True)

    no_anta_df = df[df['method'] == 'no_anta']
    methods = df['method'].unique()
    methods = [m for m in methods if m != 'no_anta']

    # Define the exclusion list
    exclusion_list = [(0.9, 0.6)]  # Add other (alpha, beta) pairs as needed

    for method in methods:
        # Matplotlib Plot
        plt.figure(figsize=(7.5, 5.5), dpi=300)
        grouped = df[df['method'] == method].groupby(['alpha', 'beta'])
        for (alpha, beta), group in grouped:
            if (alpha, beta) in exclusion_list:
                continue  # Skip the rest of the loop for this group
            no_anta_group = no_anta_df[(no_anta_df['alpha'] == alpha) & (no_anta_df['beta'] == beta)]
            sorted_group = group.sort_values('mu')
            sorted_no_anta_group = no_anta_group.sort_values('mu')
            delta_ber = sorted_group['average_ber'].values - sorted_no_anta_group['average_ber'].values
            plt.plot(sorted_group['mu'].values, delta_ber, '-o', label=f'$\\alpha$={alpha}, $\\beta$={beta}')
        # plt.title(f'$\Delta$ BER vs $\mu$ (Method: {method} - no_anta)')
        plt.xlabel(r'$\mu$', fontsize=22, labelpad=-10)
        plt.ylabel('$\Delta~BER~~~~~$', fontsize=20, labelpad=-10)
        plt.xticks([0.2, 0.3, 0.4, 0.5], fontsize=18)
        plt.yticks([0.0, -0.01, -0.02, -0.04, -0.05], fontsize=18)
        # plt.legend(loc='best')
        plt.legend(bbox_to_anchor=(1.005, 1), loc='upper left', borderaxespad=0., fontsize=14)
        plt.subplots_adjust(right=0.75)  # Adjust the bottom to make room for text annotations
        # plt.subplots_adjust(bottom=0.2)  # Adjust the bottom margin to accommodate x-axis labels

        # plt.grid(True)
        for fmt in ['png', 'tiff', 'eps']:
            matplotlib_path = os.path.join(plot_dir, f'{method}_delta_bit_error_vs_mu.{fmt}')
            plt.savefig(matplotlib_path, dpi=300, format=fmt)
        plt.close()
        # print(f'Matplotlib delta plot saved: {matplotlib_path}')

        # Plotly Plot
        fig = go.Figure()
        for (alpha, beta), group in grouped:
            no_anta_group = no_anta_df[(no_anta_df['alpha'] == alpha) & (no_anta_df['beta'] == beta)]
            sorted_group = group.sort_values('mu')
            sorted_no_anta_group = no_anta_group.sort_values('mu')
            delta_ber = sorted_group['average_ber'].values - sorted_no_anta_group['average_ber'].values
            fig.add_trace(
                go.Scatter(x=sorted_group['mu'], y=delta_ber, mode='lines+markers', name=f'Alpha={alpha}, Beta={beta}'))
        fig.update_layout(title=f'Delta BER vs $\mu$ (Method: {method} - no_anta)', xaxis_title='Mu',
                          yaxis_title='Delta BER', legend_title="Alpha, Beta")
        plotly_path = os.path.join(plot_dir, f'{method}_delta_bit_error_vs_mu.html')
        fig.write_html(plotly_path)
        # print(f'Plotly delta plot saved: {plotly_path}')


def plot_violin(csv_path):
    df = pd.read_csv(csv_path)
    plot_dir = os.path.join(os.path.dirname(csv_path), "plots")
    os.makedirs(plot_dir, exist_ok=True)

    # Filter data for anta_k_fill and no_anta methods
    filtered_df = df[df['method'].isin(['no_anta', 'anta_k_fill'])]

    # Mapping original method names to more descriptive names
    filtered_df['method'] = filtered_df['method'].map({
        'no_anta': 'No Antagonism',
        'anta_k_fill': 'Antagonism'
    })

    # Ensure that 'mu' is a categorical variable for proper ordering in the plot
    filtered_df['mu'] = pd.Categorical(filtered_df['mu'], categories=sorted(filtered_df['mu'].unique()), ordered=True)

    plt.figure(figsize=(7.5, 5.5), dpi=300)
    sns.violinplot(x='mu', y='average_ber', hue='method', data=filtered_df, split=False, palette='muted',
                   hue_order=['No Antagonism', 'Antagonism'])

    # plt.title('Violin Plot of Average BER by $\mu$')
    plt.xlabel(r'$\mu$')
    plt.ylabel('Average $BER$')
    plt.legend(loc='upper right', title='', fontsize=12)

    for fmt in ['png', 'tiff', 'eps']:
        violin_path = os.path.join(plot_dir, f'violin_plot_side_by_side_by_mu.{fmt}')
        plt.savefig(violin_path, dpi=300, format=fmt)
    plt.savefig(violin_path)
    plt.close()
    print(f'Violin plot saved: {violin_path}')


def plot_average_line(csv_path):
    df = pd.read_csv(csv_path)
    plot_dir = os.path.join(os.path.dirname(csv_path), "plots")
    os.makedirs(plot_dir, exist_ok=True)

    # Filter data for anta_k_fill and no_anta methods
    filtered_df = df[df['method'].isin(['no_anta', 'anta_k_fill'])]

    # Ensure that 'mu' is a categorical variable for proper ordering in the plot
    mu_categories = sorted(filtered_df['mu'].unique())
    filtered_df['mu'] = pd.Categorical(filtered_df['mu'], categories=mu_categories, ordered=True)

    # Calculate average bit error rate for each mu and method
    avg_ber_df = filtered_df.groupby(['mu', 'method'], as_index=False)['average_ber'].mean()

    # Separate data for each method
    no_anta_df = avg_ber_df[avg_ber_df['method'] == 'no_anta']
    anta_k_fill_df = avg_ber_df[avg_ber_df['method'] == 'anta_k_fill']

    plt.figure(figsize=(7.5, 5.5), dpi=300)

    # Plotting each method
    plt.plot(no_anta_df['mu'].values, no_anta_df['average_ber'].values, marker='o', label='no_anta')
    plt.plot(anta_k_fill_df['mu'].values, anta_k_fill_df['average_ber'].values, marker='o', label='anta_k_fill')

    plt.title('Line Plot of Average BER by $\mu$')
    plt.xlabel(r'$\mu$')
    plt.ylabel('Average BER')
    plt.legend(title='Method')

    for fmt in ['png', 'tiff', 'eps']:
        line_plot_path = os.path.join(plot_dir, f'average_line_plot_by_mu.{fmt}')
        plt.savefig(line_plot_path, dpi=300, format=fmt)
    plt.close()
    print(f'Line plot saved: {line_plot_path}')


def simulation_task(params):
    """
    Executes a single simulation task, corrects sequences using different methods,
    and saves the results along with plots.

    Parameters:
    - params: Tuple containing the parameters for the simulation (alpha, beta, mu, rho,
              num_samples, sequence_length, simulation_folder).
    """
    alpha, beta, mu, rho, num_samples, sequence_length, simulation_folder = params

    # Construct the simulation ID using rounded values
    sim_id = f"alpha_{alpha:.2f}_beta_{beta:.2f}_mu_{mu:.2f}"
    sim_path = os.path.join(simulation_folder, sim_id)
    os.makedirs(sim_path, exist_ok=True)

    print(f"Starting simulation: {sim_id}")

    # Initialize the Dataframe and apply correction methods
    df = Dataframe(alpha, beta, mu, rho, num_samples, sequence_length)
    results = []

    methods = {"viterbi": decode_with_viterbi,
               "anta_k_fill": decode_with_k,
               "no_anta": decode_with_noanta,
               "anta_k_patch": decode_with_k_patch}
    for method_name, correction_method in methods.items():
        method, average_ber = correct_and_conclude(df, correction_method, method_name, sim_path)
        results.append((alpha, beta, mu, method, average_ber))

    # Calculate the average bit error rate
    print(f"Completed: {sim_id} with Results: {results}")

    return results


def run_simulation():
    """
    Sets up and runs simulations across a range of alpha, beta, and mu values.
    Results are collected and saved to a CSV file summarizing the average BER for each set of parameters.
    """

    print("Initializing simulation parameters...")
    alphas = np.arange(0.7, 0.9, 0.1)
    betas = np.arange(0.5, 0.8, 0.1)
    mus = np.arange(0.2, 0.6, 0.1)
    # alphas = [0.7, 0.8]
    # betas = [0.5, 0.6]
    # mus = [0.2, 0.3]

    rho = 0.5
    num_samples = 1000
    sequence_length = 1000
    simulation_folder = "Simulations_Validation"
    os.makedirs(simulation_folder, exist_ok=True)

    tasks = [(round(alpha, 2), round(beta, 2), round(mu, 2), rho, num_samples, sequence_length, simulation_folder)
             for alpha in alphas for beta in betas for mu in mus]

    print("Starting simulations...")
    all_results = []
    with closing(Pool()) as pool:
        for result in pool.imap_unordered(simulation_task, tasks):
            all_results.extend(result)

    print("All simulations completed. Saving summary to CSV...")
    csv_path = os.path.join(simulation_folder, "simulation_summary.csv")
    with open(csv_path, 'w') as csv_file:
        csv_file.write("alpha,beta,mu,method,average_ber\n")
        for alpha, beta, mu, method, average_ber in all_results:
            csv_file.write(f"{alpha},{beta},{mu},{method},{average_ber}\n")

    print("Summary saved.")

    # plot_bit_error(csv_path)
    # plot_delta_error(csv_path)
    # plot_violin(csv_path)
    # plot_average_line(csv_path)


if __name__ == "__main__":
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['Arial']
    rcParams['font.size'] = 12
    # Set global font weight
    # rcParams['font.weight'] = 'bold'

    # Set global tick width
    # rcParams['axes.linewidth'] = 2  # Set the thickness of the axes lines
    rcParams['xtick.major.width'] = 2
    rcParams['ytick.major.width'] = 2
    # run_simulation()  # Comment the line if we don't want to run the simulation once again.
    # If we just want to plot it after completing the simulations
    # once to maybe change the style or so, we comment the run simulations and run the code following.
    simulation_folder = "Simulations_Validation"  # Make sure this folder exists from where the code is being run.
    csv_path = os.path.join(simulation_folder, "simulation_summary.csv")
    # plot_bit_error(csv_path)
    plot_delta_error(csv_path)
    # plot_violin(csv_path)
    # plot_average_line(csv_path)
