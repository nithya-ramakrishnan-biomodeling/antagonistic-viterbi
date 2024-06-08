import os
import pandas as pd
import multiprocessing


def process_chromosome(chromosome):
    """
    Process data for a given chromosome across various modifications.

    For each modification, this function reads the data, performs deconvolution,
    and saves the deconvolved data in a pickle file.

    Args:
    chromosome (str): The chromosome to be processed.
    """
    modifications = ['H3K27me1', 'H3K27me2', 'H3K27me3', 'H3K36me2', 'H3K36me3']
    chr_data_folder_path = 'Chr_Data'

    for modification in modifications:
        file_path = f'{chr_data_folder_path}/{chromosome}/{modification}_{chromosome}.txt'
        print(f"Processing {file_path}")

        # Read data from the file
        data = pd.read_csv(file_path, sep='\t', header=None, names=["chromosome", "start", "end", modification])
        deconvolved_list = []  # List to hold the updated rows

        # Deconvolution process
        for index, row in data.iterrows():
            start, end = row['start'], row['end']
            midpoint = (start + end) // 2

            if index == 0:
                # Special handling for the first row
                row_1 = row.copy()
                row_1['end'] = midpoint
                row_1[modification] = 0
                deconvolved_list.append(row_1)

                row_2 = row.copy()
                row_2['start'] = midpoint + 1
                row_2[modification] = 0
                deconvolved_list.append(row_2)
            else:
                # Processing subsequent rows
                row_2 = row.copy()
                row_2['start'] = midpoint + 1
                new_value = row[modification] - deconvolved_list[-1][modification]
                row_2[modification] = max(0, new_value)  # Ensure non-negative values
                deconvolved_list.append(row_2)

        # Create a DataFrame from the deconvolved list
        deconvolved_df = pd.DataFrame(deconvolved_list, columns=data.columns)
        deconvolved_df.reset_index(drop=True, inplace=True)

        # Change the file extension from .txt to .pkl for the output file
        output_file = file_path.replace('.txt', '_deconvolved.pkl')
        deconvolved_df.to_pickle(output_file)


if __name__ == '__main__':
    chromosomes = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9',
                   'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17',
                   'chr18', 'chr19', 'chrX', 'chrY', 'chrMT']

    # Use multiprocessing to process each chromosome in parallel
    with multiprocessing.Pool() as pool:
        pool.map(process_chromosome, chromosomes)
