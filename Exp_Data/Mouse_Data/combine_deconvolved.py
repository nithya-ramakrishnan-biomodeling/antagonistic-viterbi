import pandas as pd
import multiprocessing


def any_value_greater_than_0_5(row):
    """
    Check if any value in a row (starting from the fourth column) is greater than 0.5.

    Parameters:
    row (pd.Series): A row of a DataFrame.

    Returns:
    bool: True if any value in the specified part of the row is greater than 0.5, False otherwise.
    """
    return any(row[3:] > 0.5)


def mark_acceptance(row):
    """
    Determine whether a row should be marked as accepted or rejected based on specific criteria.

    Parameters:
    row (pd.Series): A row of a DataFrame.

    Returns:
    int: 1 if the row is accepted (based on the given criteria), 0 if rejected.
    """
    if any_value_greater_than_0_5(row) or (row.iloc[3:].sum() > 1.25):
        return 1  # Accepted
    else:
        return 0  # Rejected


def combine_deconvolved(chromosome):
    """
    Combine deconvolved data for a given chromosome from multiple files, apply acceptance criteria,
    and save the results.

    Parameters:
    chromosome (str): The chromosome identifier (e.g., 'chr1', 'chr2', etc.).
    """
    chr_data_folder_path = 'Chr_Data'

    # Define file paths for the specified chromosome
    file_paths = [f'{chr_data_folder_path}/{chromosome}/H3K27me1_{chromosome}_deconvolved',
                  f'{chr_data_folder_path}/{chromosome}/H3K27me2_{chromosome}_deconvolved',
                  f'{chr_data_folder_path}/{chromosome}/H3K27me3_{chromosome}_deconvolved',
                  f'{chr_data_folder_path}/{chromosome}/H3K36me2_{chromosome}_deconvolved',
                  f'{chr_data_folder_path}/{chromosome}/H3K36me3_{chromosome}_deconvolved']

    # Initialize an empty DataFrame to combine data
    combined_df = pd.DataFrame(columns=["chromosome", "start", "end"])

    # Read and merge data from each file
    for file_path in file_paths:
        print(file_path)
        file_df = pd.read_pickle(file_path + '.pkl')
        combined_df = pd.merge(combined_df, file_df, on=["chromosome", "start", "end"], how='outer')

    # Replace missing values with 'Conflict'
    combined_df = combined_df.fillna('Conflict')

    # Apply acceptance criteria
    combined_df['accepted'] = combined_df.apply(mark_acceptance, axis=1)

    # Save the data for accepted cases
    accepted_df = combined_df[combined_df['accepted'] == 1]
    accepted_df.reset_index(drop=True, inplace=True)
    accepted_df.to_pickle(f"{chr_data_folder_path}/{chromosome}/combined_accepted_data_{chromosome}.pkl")

    # Additional saving options are commented out
    # combined_df.to_pickle("combined_data_modified.pkl")
    # combined_df.to_csv("combined_data_modified.csv", index=False)
    # rejected_df = combined_df[combined_df['accepted'] == 0]
    # rejected_df.to_pickle("combined_rejected_data.pkl")
    # accepted_df.to_csv("combined_accepted_data.csv", index=False)
    # rejected_df.to_csv("combined_rejected_data.csv", index=False)

    # Print the first few rows of the accepted data (optional)
    # print(accepted_df.head())


if __name__ == '__main__':
    # Define the list of chromosomes to process
    chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
                   "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
                   "chrX", "chrY", "chrMT"]

    # Create a pool of worker processes for parallel execution
    with multiprocessing.Pool() as pool:
        # Process each chromosome in parallel
        pool.map(combine_deconvolved, chromosomes)
