import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import os
import sys
import torch
import torch.nn as nn

parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(parent_dir)

from Epigenetic_Sequence_Predictor.ESP_Sim_Viterbi import Dataframe, generate_daughter, split_sequence
from Epigenetic_Sequence_Predictor.ESP_Sim_Viterbi import viterbi_decode_antagonistic


# # Define your decoding functions here
def decode_with_viterbi(alpha, beta, rho, mu, corrupt_daughter):
    corrected_daughter = viterbi_decode_antagonistic(corrupt_daughter, alpha, beta, mu, rho)
    return corrected_daughter


def decode_with_k(alpha, beta, rho, mu, corrupt_daughter):
    corrected_sequence = []
    count = 0
    start_node = None

    k = get_k(alpha, beta, mu)

    for i in range(len(corrupt_daughter)):
        if corrupt_daughter[i] == 0:
            if count == 0:
                start_node = corrupt_daughter[i - 1] if i > 0 else None
            count += 1
        else:
            end_node = corrupt_daughter[i]
            if count > 0:
                if (start_node == 1 and end_node == 1 and count < k) or k < 0:  # if k < 0 then it is below the line.
                    # This accounts for antagonism as antagonism is a value of 2
                    corrected_sequence.extend([1] * count)
                else:
                    corrected_sequence.extend([0] * count)
                count = 0
            corrected_sequence.append(end_node)

    # Handling the case where the sequence ends with zeros
    if count > 0:
        corrected_sequence.extend([0] * count)

    # print(corrected_sequence)
    return corrected_sequence


def decode_with_k_patch(alpha, beta, rho, mu, corrupt_daughter):
    # This is the case where we use k filling with antagonism if the gap is in the vicinity of an antagonistic 1.
    # Otherwise, we use mu = 0 to calculate k and use that to fill the gap.

    corrected_sequence = []
    count = 0
    start_node = None
    vicinity_threshold = int(len(corrupt_daughter) * 0.1)  # 10 % of the length of the sequence.

    k_no_anta = get_k(alpha, beta, 0)  # k value for no antagonism
    k_anta = get_k(alpha, beta, mu)  # k value with antagonism

    for i in range(len(corrupt_daughter)):
        if corrupt_daughter[i] == 0:
            if count == 0:
                start_node = corrupt_daughter[i - 1] if i > 0 else None
                start_index = i - 1 if i > 0 else None
            count += 1
        else:
            end_node = corrupt_daughter[i]
            end_index = i
            if count > 0:
                # Adjusting the start and end indices for the vicinity region, considering boundary conditions
                vicinity_start = max(0, start_index - vicinity_threshold) if start_index is not None else 0
                vicinity_end = min(len(corrupt_daughter), end_index + vicinity_threshold)
                vicinity_region = corrupt_daughter[vicinity_start:vicinity_end]

                # Check for a 2 in the vicinity region
                if 2 not in vicinity_region:
                    # Treat is as a no antagonism case
                    k = k_no_anta
                else:
                    k = k_anta

                if (start_node == 1 and end_node == 1 and count < k) or k < 0:  # if k < 0 then it is below the line.
                    # This accounts for antagonism as antagonism is a value of 2
                    corrected_sequence.extend([1] * count)
                else:
                    corrected_sequence.extend([0] * count)
                count = 0
            corrected_sequence.append(end_node)

    # Handling the case where the sequence ends with zeros
    if count > 0:
        corrected_sequence.extend([0] * count)

    # print(corrected_sequence)
    return corrected_sequence


def decode_with_noanta(alpha, beta, rho, mu, corrupt_daughter):
    corrected_sequence = []
    count = 0
    start_node = None
    mu = 0  # Ignoring antagonism by setting mu to 0

    k = get_k(alpha, beta, mu)

    for i in range(len(corrupt_daughter)):
        if corrupt_daughter[i] == 0 or corrupt_daughter[i] == 2:  # to just look at 1s and 0s
            if count == 0:
                start_node = corrupt_daughter[i - 1] if i > 0 else None
            count += 1
        else:
            end_node = corrupt_daughter[i]
            if count > 0:
                # If the gap is between non-zero nodes (1s and 2s, or between them), fill it, ignoring k
                # As it counts it as a 0 if its 2 star_node == 1 should also work identically
                if (start_node != 0 and end_node != 0 and count < k) or k < 0:
                    corrected_sequence.extend([1] * count)
                else:
                    corrected_sequence.extend([0] * count)
                count = 0
            corrected_sequence.append(end_node)

    # Handling the case where the sequence ends with zeros
    if count > 0:
        corrected_sequence.extend([0] * count)

    # print(corrected_sequence)
    return corrected_sequence


def decode_with_nn(alpha, beta, rho, mu, corrupt_daughter, nn_model):
    # the variables alpha, beta, rho and mu are useless in our context as we are expecting a NN trained over
    # the entire range of alpha, beta and mu
    # the nn_model variable will contain the neural network with the loaded weights to work with
    corrupt_daughter_A, corrupt_daughter_B = split_sequence(corrupt_daughter)
    daughter_input = np.append(corrupt_daughter_A, corrupt_daughter_B)
    daughter_input = torch.from_numpy(daughter_input).float()

    # Run the predictions.
    nn_model.eval()  # To turn off learning mode.
    with torch.no_grad():
        corrected_daughter = nn_model(daughter_input)
    print(corrected_daughter)
    return corrected_daughter


def get_k(alpha, beta, mu):
    # Check if the point (alpha, beta) lies above the line using the new line equation
    if beta > alpha / (2 - mu):
        # Calculate k
        numerator = np.log((1 - alpha) * (1 - beta) / (alpha * beta))
        denominator = np.log((alpha / 2) / ((1 - mu / 2) * beta))
        k = np.ceil(numerator / denominator)
        # print(numerator, denominator)
        # print(f'alpha = {alpha}, beta = {beta}, mu = {mu}, k = {k}')
        return k
    else:
        # Return -99 if the point does not lie above the line
        return -99


class ExperimentalData(Dataframe):
    def __init__(self, file_paths, alpha, beta, mu, rho, mom_length):
        n_moms = len(file_paths)
        super().__init__(alpha, beta, mu, rho, n_moms, mom_length)
        self.mom_list = self._create_moms_from_files(file_paths)
        self.corrupt_daughter_list = self.mom_list.copy()
        for i, mom in enumerate(self.mom_list):
            # print(f'i = {i}')
            self.corrupt_daughter_list[i] = generate_daughter(mom, self.rho)
        self.corrected_daughter_list = self.corrupt_daughter_list.copy()

    def correct_daughter(self, method, nn_model=None):
        print(f"Correcting daughter using {method}")
        for i, corrupt_daughter in enumerate(self.corrupt_daughter_list):
            if nn_model is not None:
                corrected_daughter = method(self.alpha, self.beta, self.rho, self.mu, corrupt_daughter, nn_model)
            else:
                corrected_daughter = method(self.alpha, self.beta, self.rho, self.mu, corrupt_daughter)
            self.corrected_daughter_list[i] = corrected_daughter

    def _create_moms_from_files(self, file_paths):
        moms = []
        for file_path in file_paths:
            raw_data = np.genfromtxt(file_path, delimiter=',', skip_header=1)
            mom_sequence = self._create_mom_sequence(raw_data)
            moms.append(mom_sequence)
        return np.array(moms)

    @staticmethod
    def _create_mom_sequence(raw_data):
        mom_sequence = []
        for a_data, b_data in raw_data:
            if a_data == 0 and b_data == 0:
                mom_sequence.append(0)
            elif a_data == 1 and b_data == 0:
                mom_sequence.append(1)
            elif b_data == 1 and a_data == 0:
                mom_sequence.append(2)
            # The following condition is for using invalid sequences.
            elif a_data == 1 and b_data == 1:
                mom_sequence.append(1)
                # In this case the B is ignored.
            else:
                raise ValueError('Invalid data encountered in sequence')
        return np.array(mom_sequence)


class SequencePredictor(nn.Module):
    # Define the neural network architecture
    # (input -> h1 -> h2 .... hn -> output)
    def __init__(self):
        super(SequencePredictor, self).__init__()
        self.input_size = 200
        self.output_size = 100
        self.hidden_size = 150
        self.num_hid_layers = 5

        # This is the first layer which takes in inputs and gives out the outputs for the first hidden layer
        self.input_layer = nn.Linear(self.input_size, self.hidden_size)

        # These are the hidden layers
        self.hidden_layers = nn.ModuleList(
            [nn.Linear(self.hidden_size, self.hidden_size) for _ in range(self.num_hid_layers)])

        # This is the output layer
        self.output_layer = nn.Linear(self.hidden_size, self.output_size)

        self.relu = nn.ReLU()

    def forward(self, x):
        out = self.input_layer(x)
        out = self.relu(out)
        for layer in self.hidden_layers:
            out = layer(out)
            out = self.relu(out)

        out = self.output_layer(out)
        return out


def plot_sequences(sim_path,
                   plot_mom_a=True, plot_mom_b=True,
                   plot_corrupt_a=True, plot_corrupt_b=True,
                   plot_corrected_a=True, plot_corrected_b=True,
                   start_index=0, end_index=None):
    # Paths to the files generated by conclusions function
    mom_file = os.path.join(sim_path, 'Mom_list.csv')
    corrupt_daughter_file = os.path.join(sim_path, 'Corrupt_Daughter.csv')
    correct_daughter_file = os.path.join(sim_path, 'Correct_Daughter.csv')

    # Load the data
    mom_data = np.loadtxt(mom_file, delimiter=',')
    corrupt_daughter_data = np.loadtxt(corrupt_daughter_file, delimiter=',')
    correct_daughter_data = np.loadtxt(correct_daughter_file, delimiter=',')

    # Split sequences into A and B components
    mom_data_A, mom_data_B = split_sequence(mom_data[0])
    corrupt_daughter_data_A, corrupt_daughter_data_B = split_sequence(corrupt_daughter_data[0])
    correct_daughter_data_A, correct_daughter_data_B = split_sequence(correct_daughter_data[0])

    # If end_index is not provided, plot till the end of the sequence
    if end_index is None:
        end_index = mom_data.shape[1]

    # Ensure valid index range
    start_index = max(0, start_index)
    end_index = min(mom_data.shape[1], end_index)

    # Determine which plots to include
    plots_to_include = [
        plot_mom_a, plot_mom_b,
        plot_corrupt_a, plot_corrupt_b,
        plot_corrected_a, plot_corrected_b
    ]

    data_sets = [
        mom_data_A, mom_data_B,
        corrupt_daughter_data_A, corrupt_daughter_data_B,
        correct_daughter_data_A, correct_daughter_data_B
    ]

    titles = ['Mom A', 'Mom B', 'Corrupt A', 'Corrupt B', 'Corrected A', 'Corrected B']
    colors = ['blue', 'red', 'blue', 'red', 'blue', 'red']

    # Plotting only the specified ranges for the first sequence
    fig, axes = plt.subplots(nrows=sum(plots_to_include), ncols=1, figsize=(12, 3 * sum(plots_to_include)), sharex=True)

    # Create Plotly subplot figure
    plotly_fig = make_subplots(rows=sum(plots_to_include), cols=1, shared_xaxes=True)

    ax_idx = 0
    for i, (include, data, title, color) in enumerate(zip(plots_to_include, data_sets, titles, colors)):
        if include:
            axes[ax_idx].bar(range(start_index, end_index), data[start_index:end_index], color=color)
            axes[ax_idx].set_ylabel('Value')
            axes[ax_idx].text(1.02, 0.5, title, va='center', ha='left', rotation='vertical',
                              transform=axes[ax_idx].transAxes)
            ax_idx += 1

    if ax_idx > 0:
        axes[-1].set_xlabel('Index')
        plt.subplots_adjust(hspace=0.3)

        # Save the Matplotlib figure
        plt.savefig(f"{sim_path}/plot.png")
        print(f'Saved at {sim_path}/plot.png')

        # plt.show()

        # Close the figure to release memory
        plt.close(fig)

    # Plotly plotting code
    row_idx = 1
    for i, (include, data, title, color) in enumerate(zip(plots_to_include, data_sets, titles, colors)):
        if include:
            trace = go.Bar(x=list(range(start_index, end_index)), y=data[start_index:end_index],
                           name=title, marker_color=color)
            plotly_fig.add_trace(trace, row=row_idx, col=1)
            row_idx += 1

    # Update layout and save the Plotly figure
    plotly_fig.update_layout(height=300 * sum(plots_to_include), title_text="Sequence Data")
    plotly_fig.write_html(f"{sim_path}/plot.html")


def validate_sequence(chromosome, filter_window, binarize_threshold, sequence_number, num_samples=10, mom_length=1000,
                      chr_data_folder_path='Chr_Data', network_weight_path='extended_network_weights.pth', rho=0.5,
                      make_plots=True,
                      plot_mom_a=True, plot_mom_b=True,
                      plot_corrupt_a=True, plot_corrupt_b=False,
                      plot_corrected_a=True, plot_corrected_b=False,
                      start_index=0, end_index=None):
    # print('Starting Validation.')
    sequence_detection_path = os.path.join(chr_data_folder_path, chromosome,
                                           f"sequence_detection_{filter_window}_{binarize_threshold}")

    all_sequences_file = os.path.join(sequence_detection_path, "all_sequences.csv")

    sequence_file = os.path.join(sequence_detection_path, f"all_sequences/sequence_{sequence_number}.csv")

    # Read all_sequence.csv into a DataFrame
    all_sequences_df = pd.read_csv(all_sequences_file)
    # print(all_sequences_df)
    # print('sequence_number = ', sequence_number)
    # print(all_sequences_df['Sequence_number'] == sequence_number)
    # Find the row in the DataFrame corresponding to the sequence number
    sequence_row = all_sequences_df[all_sequences_df['Sequence_number'] == sequence_number]
    # Find alpha, beta and mu values for the sequence
    # print('sequence_row = ', sequence_row)

    alpha = sequence_row['Alpha'].iloc[0]
    beta = sequence_row['Beta'].iloc[0]
    mu = sequence_row['Mu'].iloc[0]

    n_sequence_file = [sequence_file for _ in range(num_samples)]
    exp_data = ExperimentalData(n_sequence_file, alpha, beta, mu, rho, mom_length)

    # Solve using Viterbi and store results
    viterbi_save_path = os.path.join(sequence_detection_path, f"all_sequences/sequence_{sequence_number}",
                                     "viterbi_decode")

    if not os.path.exists(viterbi_save_path):
        os.makedirs(viterbi_save_path, exist_ok=True)
        exp_data.correct_daughter(decode_with_viterbi)
        exp_data.conclusions(viterbi_save_path, level=1)
        if make_plots:
            plot_sequences(viterbi_save_path, plot_mom_a, plot_mom_b, plot_corrupt_a, plot_corrupt_b, plot_corrected_a,
                           plot_corrected_b, start_index, end_index)
    else:
        # if folder exists, check for plots, if no plots, analysis is already done, so create the plots.
        # we only have to check one plot. Assume the other is there.
        plot_path = os.path.join(viterbi_save_path, 'plot.html')
        if not os.path.exists(plot_path):
            if make_plots:
                plot_sequences(viterbi_save_path, plot_mom_a, plot_mom_b, plot_corrupt_a, plot_corrupt_b,
                               plot_corrected_a,
                               plot_corrected_b, start_index, end_index)

    # Solve using K-Fill and store results
    k_fill_save_path = os.path.join(sequence_detection_path, f"all_sequences/sequence_{sequence_number}",
                                    "anta_k_fill_decode")

    if not os.path.exists(k_fill_save_path):
        os.makedirs(k_fill_save_path, exist_ok=True)
        exp_data.correct_daughter(decode_with_k)
        exp_data.conclusions(k_fill_save_path, level=1)
        if make_plots:
            plot_sequences(k_fill_save_path, plot_mom_a, plot_mom_b, plot_corrupt_a, plot_corrupt_b, plot_corrected_a,
                           plot_corrected_b, start_index, end_index)
    else:
        # if folder exists, check for plots, if no plots, analysis is already done, so create the plots.
        # we only have to check one plot. Assume the other is there.
        plot_path = os.path.join(k_fill_save_path, 'plot.html')
        if not os.path.exists(plot_path):
            if make_plots:
                plot_sequences(k_fill_save_path, plot_mom_a, plot_mom_b, plot_corrupt_a, plot_corrupt_b,
                               plot_corrected_a,
                               plot_corrected_b, start_index, end_index)

    # No Antagonism Decoding
    no_anta_save_path = os.path.join(sequence_detection_path, f"all_sequences/sequence_{sequence_number}",
                                     "no_anta_decode")
    if not os.path.exists(no_anta_save_path):
        os.makedirs(no_anta_save_path, exist_ok=True)
        exp_data.correct_daughter(decode_with_noanta)
        exp_data.conclusions(no_anta_save_path, level=1)
        if make_plots:
            plot_sequences(no_anta_save_path, plot_mom_a, plot_mom_b, plot_corrupt_a, plot_corrupt_b, plot_corrected_a,
                           plot_corrected_b, start_index, end_index)
    else:
        # if folder exists, check for plots, if no plots, analysis is already done, so create the plots.
        # we only have to check one plot. Assume the other is there.
        plot_path = os.path.join(no_anta_save_path, 'plot.html')
        if not os.path.exists(plot_path):
            if make_plots:
                plot_sequences(no_anta_save_path, plot_mom_a, plot_mom_b, plot_corrupt_a, plot_corrupt_b,
                               plot_corrected_a,
                               plot_corrected_b, start_index, end_index)

    # Solve using K-Fill with Patch and store results
    k_patch_save_path = os.path.join(sequence_detection_path, f"all_sequences/sequence_{sequence_number}",
                                     "anta_k_patch_decode")

    if not os.path.exists(k_patch_save_path):
        os.makedirs(k_patch_save_path, exist_ok=True)
        exp_data.correct_daughter(decode_with_k_patch)
        exp_data.conclusions(k_patch_save_path, level=1)
        if make_plots:
            plot_sequences(k_patch_save_path, plot_mom_a, plot_mom_b, plot_corrupt_a, plot_corrupt_b, plot_corrected_a,
                           plot_corrected_b, start_index, end_index)
    else:
        # if folder exists, check for plots, if no plots, analysis is already done, so create the plots.
        # we only have to check one plot. Assume the other is there.
        plot_path = os.path.join(k_patch_save_path, 'plot.html')
        if not os.path.exists(plot_path):
            if make_plots:
                plot_sequences(k_patch_save_path, plot_mom_a, plot_mom_b, plot_corrupt_a, plot_corrupt_b,
                               plot_corrected_a,
                               plot_corrected_b, start_index, end_index)

    # Solve using NN and store results
    nn_save_path = os.path.join(sequence_detection_path, f"all_sequences/sequence_{sequence_number}",
                                "nn_decode")
    model = SequencePredictor()
    model.load_state_dict(torch.load(network_weight_path))  # Load the network weights

    if not os.path.exists(nn_save_path):
        os.makedirs(nn_save_path, exist_ok=True)
        exp_data.correct_daughter(decode_with_nn, model)
        exp_data.conclusions(nn_save_path, level=1)
        if make_plots:
            plot_sequences(nn_save_path, plot_mom_a, plot_mom_b, plot_corrupt_a, plot_corrupt_b, plot_corrected_a,
                           plot_corrected_b, start_index, end_index)
        else:
            # if folder exists, check for plots, if no plots, analysis is already done, so create the plots.
            # we only have to check one plot. Assume the other is there.
            plot_path = os.path.join(nn_save_path, 'plot.html')
            if not os.path.exists(plot_path):
                if make_plots:
                    plot_sequences(nn_save_path, plot_mom_a, plot_mom_b, plot_corrupt_a, plot_corrupt_b,
                                   plot_corrected_a,
                                   plot_corrected_b, start_index, end_index)


if __name__ == "__main__":
    validate_sequence(chromosome='chr1', filter_window='5', binarize_threshold='05', sequence_number=170539,
                      num_samples=10,
                      mom_length=1000,
                      make_plots=True,
                      plot_mom_a=True, plot_mom_b=True,
                      plot_corrupt_a=True, plot_corrupt_b=False,
                      plot_corrected_a=True, plot_corrected_b=False,
                      start_index=0, end_index=None)


# Supplementary

def plot_line_and_point(alpha, beta, mu):
    # Create a range of alpha values
    alpha_range = np.linspace(0, 1, 500)

    # Calculate beta values for the line using the new equation
    beta_line = alpha_range / (2 - mu)

    # Plot the line
    plt.plot(alpha_range, beta_line, label='Boundary Line')

    # Plot the point
    plt.scatter([alpha], [beta], color='red')
    plt.text(alpha, beta, f'  ({alpha}, {beta})')

    # Setting plot limits
    plt.xlim(0, 1)
    plt.ylim(0, 1)

    # Adding labels and title
    plt.xlabel('Alpha')
    plt.ylabel('Beta')
    plt.title('Plot of Line and Point')
    plt.legend()

    # Show the plot
    plt.show()
