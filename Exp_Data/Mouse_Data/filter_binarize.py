import pandas as pd
# from scipy.signal import medfilt
import numpy as np
import multiprocessing
import os


def binarize_data(data, columns_to_binarize, threshold=0.5):
    """
    Binarize the data based on a given threshold.
    Args:
    data (DataFrame): The pandas DataFrame containing the data.
    columns_to_binarize (list): List of column indices to binarize.
    threshold (float): The threshold for binarization.
    Returns:
    DataFrame: The binarized data.
    """
    for col_index in columns_to_binarize:
        # Use apply with lambda for binarization and explicit conversion to int
        data.iloc[:, col_index] = data.iloc[:, col_index].apply(lambda x: 1 if x > threshold else 0).astype(int)
        col_name = data.columns[col_index]
        data[col_name] = data[col_name].astype(int)
    return data


def median_filter_data(data, columns_to_filter, window_size=5):
    """
    Apply median filtering to the data.
    Args:
    data (DataFrame): The pandas DataFrame containing the data.
    columns_to_filter (list): List of column indices to apply median filtering.
    window_size (int): The size of the median filter window.
    Returns:
    DataFrame: The filtered data.
    """
    # In case Scipy doesnt work
    # for col in columns_to_filter:
    #     data.iloc[:, col] = medfilt(data.iloc[:, col], kernel_size=window_size)

    half_window = window_size // 2

    for col in columns_to_filter:
        col_data = data.iloc[:, col].tolist()
        filtered_col = []

        for i in range(len(col_data)):
            # Calculate window boundaries
            start_index = max(i - half_window, 0)
            end_index = min(i + half_window + 1, len(col_data))

            # Extract the window and compute the median
            window = col_data[start_index:end_index]
            median_value = np.median(window)
            filtered_col.append(median_value)

        # Replace the original column data with the filtered data
        data.iloc[:, col] = filtered_col

    return data


def process_data(chromosome):
    """
    Process data for a given chromosome.
    Perform median filtering and then binarize the data, saving the result.
    Args:
    chromosome (str): The chromosome to be processed.
    """
    chr_data_folder_path = 'Chr_Data'
    window_size = 5
    binarize_threshold = 0.5

    binarize_threshold_str = str(binarize_threshold).replace('.', '')

    # Read the combined file
    input_file = f"{chr_data_folder_path}/{chromosome}/combined_accepted_data_{chromosome}.pkl"
    print(input_file)

    combined_data = pd.read_pickle(input_file)

    # Define column numbers to modify
    columns_to_modify = [3, 4, 5, 6, 7, 8]  # Replace with appropriate column indices

    # Perform median filtering and then binarize data
    filtered_data = median_filter_data(combined_data.copy(), columns_to_modify, window_size)
    filtered_then_binarized_data = binarize_data(filtered_data, columns_to_modify, binarize_threshold)
    filtered_then_binarized_data.reset_index(drop=True, inplace=True)
    filtered_then_binarized_data.to_pickle(
        f"{chr_data_folder_path}/{chromosome}/filtered_{window_size}_then_binarized_{binarize_threshold_str}_{chromosome}.pkl")

    # Other processing options (commented out for potential future use)
    # binarized_data = binarize_data(combined_data.copy(), columns_to_modify, binarize_threshold)
    # binarized_then_filtered_data = median_filter_data(binarized_data, columns_to_modify, window_size)
    # binarized_then_filtered_data.reset_index(drop=True, inplace=True)
    # binarized_then_filtered_data.to_pickle(
    # f"{chr_data_folder_path}/{chromosome}/binarized_{binarize_threshold_str}_then_filtered_{window_size}_{chromosome}.pkl")


if __name__ == '__main__':
    chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
                   "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
                   "chrX", "chrY", "chrMT"]

    with multiprocessing.Pool() as pool:
        pool.map(process_data, chromosomes)
