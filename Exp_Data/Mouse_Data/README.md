# Experimental Validation

## Working with the Experimental Data

* The experimental data collected from the Groth-Lab Papers are utilized for experimental validation.
* To investigate antagonistic traits, special attention is given to the data on H3K27me3 and H3K36me3.
* As the data do not directly test for the presence of H3 alone, it is necessary to deduce which nucleosomes contain H3.
* This is accomplished by considering other modifications as well and then employing a sorting algorithm to identify
  nucleosomes with H3.

## Data Used

* GSM3983930_H3K27me1_ChIP-seq_rep1.txt -> Renamed to H3K27me1.txt
* GSM3983932_H3K27me2_ChIP-seq_rep1.txt -> Renamed to H3K27me2.txt
* GSM3983934_H3K27me3_ChIP-seq_rep1.txt -> Renamed to H3K27me3.txt
* GSM3983936_H3K36me2_ChIP-seq_rep1.txt -> Renamed to H3K36me2.txt
* GSM3983938_H3K36me3_ChIP-seq_rep1.txt -> Renamed to H3K36me3.txt

These files have been renamed for simplicity and are stored in the `Raw_Data` directory within this workspace.

## Chromosome Based Sorting

* The raw data encompasses chromosomes {1 - 19, MT, X, and Y}, resulting in a considerable size that is challenging to
  manage.
* To address this, the raw data is sorted into individual folders, one for each chromosome.
* The sorting is executed using the `sort_by_chr.py` Python script.
* This script inputs data from the `Raw_Data` folder and outputs a new folder named `Chr_Data`.
* This approach allows for more manageable work on individual chromosomes, which are smaller in scale compared to the
  complete dataset.

## Processing Chromosomal Data

With the data now segmented into smaller chromosome-specific portions, the following analyses can be conducted on each
chromosome. Hence, references to "data" in the ensuing discussion pertain to chromosomal data.

The majority of program outputs are in `.pkl` format to expedite calculations, though options exist to save in `.csv`
format if needed (though not recommended). Use the `pandas` library for converting `.pkl` files to `.csv`, if necessary.

## Deconvolution

* Each modification in the data specifies a start and end base pair number. Notably, there is an overlap of 500 base
  pairs in two consecutive rows.
* This overlap necessitates the use of a deconvolution algorithm to extract the required nucleosome data.
* The `deconvolve_raw` Python script performs this task, processing data in the `Chr_Data` path.
* It processes each chromosome folder, converting raw files into their deconvolved forms.
* On servers with Slurm, use the `run_deconvolve.sh` script with the `sbatch` command.
* Ensure the `.py`, `.sh` scripts, and the `Chr_Data` folder are in the same directory and execute the `.sh` script from
  there.

## Combine and Sort

* In the deconvolved file, each row (500 base pairs) is considered a nucleosome. Practically, a nucleosome is about 200
  base pairs, but for error minimization, we assume two nucleosomes combined act as one.
* The presence of H3 is determined using a combination of all data files and a sorting criterion.
* The sorting criteria are as follows: H3 is presumed present if
    1. Any modification's value exceeds 0.5 or
    2. The sum of all modification values is greater than 1.25.
* This process is handled by the `combine_deconvolved.py` Python script.
* The script processes each chromosome in the `Chr_Data` folder, generating combined and sorted data according to the
  above criteria.
* Additional script functionalities are commented out but available for deeper data analysis.
* The output is a `.pkl` file for each chromosome, containing combined and sorted data. Other selection criteria can be
  defined in the script.
* For Slurm servers, use the `run_combine.sh` script with `sbatch`.
* Place the `.py`, `.sh` scripts, and the `Chr_Data` folder in the same directory and run the `.sh` script from there.

## Filter and Binarize

* Our study focuses on H3K27me3 and H3K36me3, known antagonistic modifications, in nucleosomes containing H3.
* The data often exhibits abrupt value changes, requiring a median filter (window size of 5) for smoothing.
* Post-filtering, the data is binarized using a 0.5 threshold. Alternatively, binarize first, then apply the median
  filter.
* The `filter_binarize.py` Python script performs these tasks.
* It processes the `Chr_Data` folder, filtering and binarizing the combined and sorted data for all modifications.
* Additional script functionalities include options for first binarizing, then filtering, or providing only filtered
  data.
* The output is a `.pkl` file for each chromosome, labeled with filter window size and binarization threshold details.
* For instance, `filtered_5_then_binarized_05_chr1.pkl` indicates chr1 data filtered with a window of 5 and binarized at
  0.5.
* For Slurm servers, use the `run_filter.sh` script with `sbatch`.
* Place the `.py`, `.sh` scripts, and the `Chr_Data` folder in the same directory and execute the `.sh` script from
  there. Here's an enhanced and clarified version of the "Sequence Detection" section of your documentation:

## Sequence Detection - All

* The binarized dataset is extensive, containing millions of nucleosomes. For algorithmic verification, it is crucial to
  work with smaller sequences.
* We aim to select sequence subsets from the entire dataset for validation purposes.
* Our theoretical framework is applicable to specific regions within the 3D space defined by alpha, beta, and mu
  parameters.
* We first create all possible sequences and id them such that they can be picked out based on conditions we can later
  decide.
* To identify optimal regions of interest, sequences of 1000 nucleosomes are scanned with overlapping segments.
* The `detect_seq.py` Python script facilitates this process.
* It requires the `Chr_Data` path and the filter window and binarization threshold values to accurately select the file
  generated by the `filter_binarize` program.
* The script processes each chromosome's folder, utilizing the filtered and binarized data to compile a list of all
  sequences.
* Each row in the output list corresponds to one sequence.
* A folder named `sequence_detection_{filter_window}_{binarization_threshold}` is generated, where all data specific to
  that filter window and binarization threshold are stored. For example, with a filter window of 5 and a binarization
  threshold of 0.5, the folder would be `sequence_detection_5_05`.
* All subsequent files and folders are organized within this directory.
* The output file `all_sequences.csv`, for each chromosome, contains a list of all sequences from that specific dataset.
  This file remains consistent across all selection criteria for a given dataset.
* A `all_sequences` folder is created to store all sequences, labeled by their sequence number.
* Detailed information about each selected sequence, like its starting base pair, is provided in the `all_sequences.csv`
  file, associated with its sequence number.
* For execution on servers with Slurm, utilize the `run_detect_all.sh` script with the `sbatch` command.
* Ensure that the `.py`, `.sh` scripts, and the `Chr_Data` folder are all located in the same directory, and execute
  the `.sh` script from that directory.

## Data Validation

* Having created all the sequences, we need functions that can take in the sequence, validate it using different
  methods, like viterbi, k-fill etc..
* This is done using the `data_validation.py` program.
* 

## Sequence Selection - Webserver

* Because the number of files in the Chr_Data file is too large, it is easier to host a webserver to show us all the
  files and lets us filter through the parameters that we want. It will also allow us to validate the data that we
  select such that
* The webserver and all necessary files are stored in the `sequence_webserver` directory

* A flask program that creates a server that takes in the `all_sequences.csv` and creates a website that lets us sort
  and select sequences from all that has been created.