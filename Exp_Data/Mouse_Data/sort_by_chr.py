import pandas as pd
import os
import multiprocessing


def process_modification(mod):
    """
    Process the data for a given modification type.

    This function reads the data from a file corresponding to the modification type,
    groups the data by chromosome, and then saves each group in a separate file
    in a directory named after the chromosome.

    Args:
    mod (str): The modification type to process.
    """
    # Define the data folder paths within the function
    raw_data_folder_path = 'Raw_Data'
    chr_data_folder_path = 'Chr_Data'

    file_name = f'{raw_data_folder_path}/{mod}.txt'
    print(f"Processing {file_name}")

    column_names = ['Chromosome', 'Start', 'End', 'Value']
    data = pd.read_csv(file_name, sep='\t', header=None, names=column_names)

    # Group data by Chromosome
    grouped_data = data.groupby('Chromosome')

    for chromosome, chromosome_data in grouped_data:
        # Create directory for each chromosome within Chr_Data if it doesn't exist
        chromosome_folder = os.path.join(chr_data_folder_path, chromosome)
        if not os.path.exists(chromosome_folder):
            os.makedirs(chromosome_folder)

        # Save the grouped data in its corresponding chromosome folder
        output_file = os.path.join(chromosome_folder, f"{mod}_{chromosome}.txt")
        chromosome_data.to_csv(output_file, sep='\t', header=False, index=False)


if __name__ == '__main__':
    modifications = ['H3K27me1', 'H3K27me2', 'H3K27me3', 'H3K36me2', 'H3K36me3']

    # Use multiprocessing to process each modification in parallel
    with multiprocessing.Pool() as pool:
        pool.map(process_modification, modifications)
