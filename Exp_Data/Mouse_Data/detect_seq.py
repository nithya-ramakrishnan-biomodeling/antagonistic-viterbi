import pandas as pd
import multiprocessing
import os


def calculate_alpha_beta(binary_sequence):
    """
    Calculates alpha and beta values for a given binary sequence.
    Alpha is the transition probability from 1 to 1, and Beta is from 0 to 0.

    Parameters:
    binary_sequence (str): A string of binary digits (1s and 0s).

    Returns:
    tuple: A tuple containing the alpha and beta values.
    """
    n11, n10, n01, n00 = 0, 0, 0, 0
    sequence = [int(bit) for bit in binary_sequence]

    for i in range(1, len(sequence)):
        if sequence[i - 1] == 1 and sequence[i] == 1:
            n11 += 1
        elif sequence[i - 1] == 1 and sequence[i] == 0:
            n10 += 1
        elif sequence[i - 1] == 0 and sequence[i] == 1:
            n01 += 1
        elif sequence[i - 1] == 0 and sequence[i] == 0:
            n00 += 1

    alpha = n11 / (n11 + n10) if (n11 + n10) > 0 else 0
    beta = n00 / (n00 + n01) if (n00 + n01) > 0 else 0

    return alpha, beta


def calculate_mu(a_data, b_data):
    """
    Calculates the probability Mu, the probability that B_data is 1 given A_data is 0.

    Parameters:
    a_data (list): A list of binary values representing A_data.
    b_data (list): A list of binary values representing B_data.

    Returns:
    float: The calculated probability Mu.
    """
    count_a0_b1 = sum(1 for a, b in zip(a_data, b_data) if a == 0 and b == 1)
    count_a0 = a_data.count(0)

    mu = count_a0_b1 / count_a0 if count_a0 > 0 else 0
    return mu


def detect_sequence(chromosome):
    """
    Processes each chromosome to detect and filter sequences based on specific criteria.
    The function reads filtered and binarized data and saves valid sequences and a summary of results.

    Parameters:
    chromosome (str): The chromosome to process.
    """
    print(chromosome)
    chr_data_folder_path = 'Chr_Data'
    filter_window = 5
    binarize_threshold = 0.5

    binarize_threshold_str = str(binarize_threshold).replace('.', '')
    file_name = f"{chr_data_folder_path}/{chromosome}/filtered_{filter_window}_then_binarized_{binarize_threshold_str}_{chromosome}.pkl"
    df = pd.read_pickle(file_name)

    a_data = df['H3K27me3']
    b_data = df['H3K36me3']
    start_positions = df['start']
    end_positions = df['end']

    # New directory structure
    sequence_detection_dir = f"{chr_data_folder_path}/{chromosome}/sequence_detection_{filter_window}_{binarize_threshold_str}"
    all_sequences_dir = f"{sequence_detection_dir}/all_sequences"
    os.makedirs(all_sequences_dir, exist_ok=True)

    all_sequences = []

    batch_size = 1000

    for i in range(len(a_data) - (batch_size - 1)):
        if i % 10000 == 0 or i == len(a_data) - 10000:
            print(f'{chromosome} :: {i}/{len(a_data) - 999}')

        window_a = a_data[i:i + batch_size]
        window_b = b_data[i:i + batch_size]
        start_bp = start_positions[i]
        end_bp = end_positions[i]

        mu = calculate_mu(window_a.tolist(), window_b.tolist())
        alpha, beta = calculate_alpha_beta(''.join(map(str, window_a)))
        invalid = sum(1 for a, b in zip(window_a, window_b) if a == 1 and b == 1)

        file_name = f"{all_sequences_dir}/sequence_{i}.csv"
        pd.DataFrame({'A_data': window_a, 'B_data': window_b}).to_csv(file_name, index=False)

        all_sequences.append([i, start_bp, end_bp, alpha, beta, mu, invalid])

    all_sequences_df = pd.DataFrame(all_sequences,
                                    columns=['Sequence_number', 'Start_BP', 'End_BP', 'Alpha', 'Beta', 'Mu', 'Invalid'])
    all_sequences_df.to_csv(f"{sequence_detection_dir}/all_sequences.csv", index=False)


if __name__ == '__main__':
    chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
                   "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
                   "chrX", "chrY", "chrMT"]

    with multiprocessing.Pool() as pool:
        pool.map(detect_sequence, chromosomes)

    print('Completed!')
