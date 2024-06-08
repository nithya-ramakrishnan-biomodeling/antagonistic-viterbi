from flask import Flask, render_template, request, send_file
import pandas as pd
import os
import sys
from multiprocessing import Pool

parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(parent_dir)

from validate_data import validate_sequence, get_k

app = Flask(__name__)

CHR_DATA_PATH = "../Chr_Data"
NN_WEIGHTS = "../extended_network_weights.pth"


@app.route('/', methods=['GET', 'POST'])
def index():
    chromosomes = [_ for _ in os.listdir(CHR_DATA_PATH) if os.path.isdir(os.path.join(CHR_DATA_PATH, _))]

    # Assuming the first chromosome's folders are representative of all
    first_chr_path = os.path.join(CHR_DATA_PATH, "chr1")
    sequence_detection_folders = [folder for folder in os.listdir(first_chr_path) if "sequence_detection" in folder]
    windows = set(folder.split('_')[2] for folder in sequence_detection_folders)
    thresholds = set(folder.split('_')[3] for folder in sequence_detection_folders)

    # Initialize variables
    filtered_data = main_data = None
    total = 0
    current_page = 1
    rows_per_page = 20

    # Handle form submission
    if request.method == 'POST':
        chr_selected = request.form.get('chromosome')
        window_selected = request.form.get('window')
        threshold_selected = request.form.get('threshold')
        alpha_min = float(request.form.get('alpha_min', 0))
        alpha_max = float(request.form.get('alpha_max', 1))
        beta_min = float(request.form.get('beta_min', 0))
        beta_max = float(request.form.get('beta_max', 1))
        mu_min = float(request.form.get('mu_min', 0))
        mu_max = float(request.form.get('mu_max', 1))
        invalid_min = int(request.form.get('invalid_min', 0))
        invalid_max = int(request.form.get('invalid_max', 100))
        current_page = int(request.form.get('current_page', 1))

        folder_name = f"sequence_detection_{window_selected}_{threshold_selected}"
        file_path = os.path.join(CHR_DATA_PATH, chr_selected, folder_name, 'all_sequences.csv')

        if os.path.exists(file_path):
            df = pd.read_csv(file_path)
            filtered_data = df[(df['Alpha'] >= alpha_min) & (df['Alpha'] <= alpha_max) &
                               (df['Beta'] >= beta_min) & (df['Beta'] <= beta_max) &
                               (df['Mu'] >= mu_min) & (df['Mu'] <= mu_max) &
                               (df['Invalid'] >= invalid_min) & (df['Invalid'] <= invalid_max)]
            total = len(filtered_data)

            # Pagination logic
            start = (current_page - 1) * rows_per_page
            end = start + rows_per_page
            main_data = filtered_data.iloc[start:end]

    return render_template('index.html',
                           chromosomes=chromosomes,
                           windows=windows,
                           thresholds=thresholds,
                           data=main_data,
                           total=total,
                           current_page=current_page,
                           rows_per_page=rows_per_page)


def perform_analysis(args, chr_data_path=CHR_DATA_PATH, network_weight_path=NN_WEIGHTS):
    index, chromosome, window, threshold, sequence_number, make_plots = args
    print(f'Working on {index} index')
    folder_name = f"sequence_detection_{window}_{threshold}"
    file_path = os.path.join(chr_data_path, chromosome, folder_name, 'all_sequences.csv')
    if not os.path.exists(file_path):
        return None, "File not found"

    df = pd.read_csv(file_path)
    sequence_data = df[df['Sequence_number'] == int(sequence_number)].iloc[0]
    validate_sequence(chromosome=chromosome, filter_window=window, binarize_threshold=threshold,
                      sequence_number=int(sequence_number), num_samples=10, mom_length=200,
                      chr_data_folder_path=chr_data_path, network_weight_path=network_weight_path, rho=0.5,
                      make_plots=make_plots, plot_mom_a=True,
                      plot_mom_b=True, plot_corrupt_a=True, plot_corrupt_b=True, plot_corrected_a=True,
                      plot_corrected_b=False, start_index=0, end_index=None)

    decode_types = ['anta_k_fill', 'anta_k_patch', 'viterbi', 'no_anta', 'nn']
    avg_bit_errors = {}
    bit_error_deltas = {}

    for decode_type in decode_types:
        bit_error_file_path = os.path.join(chr_data_path, chromosome, folder_name, 'all_sequences',
                                           f"sequence_{int(sequence_number)}", f"{decode_type}_decode", "BitError.csv")
        if os.path.exists(bit_error_file_path):
            bit_error_df = pd.read_csv(bit_error_file_path, header=None)
            avg_bit_errors[decode_type] = bit_error_df.iloc[:, 0].mean()
        else:
            avg_bit_errors[decode_type] = "N/A"

    no_anta_avg = avg_bit_errors.get('no_anta', 0)

    for decode_type in decode_types:
        if decode_type != 'no_anta' and avg_bit_errors[decode_type] != "N/A":
            bit_error_deltas[decode_type] = avg_bit_errors[decode_type] - no_anta_avg

    k_no_antagonism = get_k(sequence_data['Alpha'], sequence_data['Beta'], 0)
    k_antagonism = get_k(sequence_data['Alpha'], sequence_data['Beta'], sequence_data['Mu'])

    sequence_details = {
        'Sequence_Number': sequence_number,
        'Start_BP': sequence_data['Start_BP'],
        'End_BP': sequence_data['End_BP'],
        'Alpha': sequence_data['Alpha'],
        'Beta': sequence_data['Beta'],
        'Mu': sequence_data['Mu'],
        'Invalid': sequence_data['Invalid'],
        'K_Value_Antagonism': k_antagonism,
        'K_Value_No_Antagonism': k_no_antagonism,
        'Avg_Bit_Error_K_Fill_Antagonism': avg_bit_errors.get('anta_k_fill', 'N/A'),
        'Avg_Bit_Error_K_Patch_Antagonism': avg_bit_errors.get('anta_k_patch', 'N/A'),
        'Avg_Bit_Error_Viterbi': avg_bit_errors.get('viterbi', 'N/A'),
        'Avg_Bit_Error_No_Anta': avg_bit_errors.get('no_anta', 'N/A'),
        'Avg_Bit_Error_NN': avg_bit_errors.get('nn', 'N/A'),
        'Delta_Bit_Error_K_Fill_Antagonism': bit_error_deltas.get('anta_k_fill', 0),
        'Delta_Bit_Error_K_Patch_Antagonism': bit_error_deltas.get('anta_k_patch', 0),
        'Delta_Bit_Error_Viterbi': bit_error_deltas.get('viterbi', 0),
        'Delta_Bit_Error_NN': bit_error_deltas.get('nn', 0)
    }

    return {'file_path': file_path,
            'chromosome': chromosome,
            'window': window,
            'threshold': threshold,
            'sequence_data': sequence_data,
            'avg_bit_errors': avg_bit_errors,
            'bit_error_deltas': bit_error_deltas,
            'k_no_antagonism': k_no_antagonism,
            'k_antagonism': k_antagonism,
            }, sequence_details, None


@app.route('/analyze', methods=['POST'])
def analyze():
    chromosome = request.form.get('chromosome')
    window = request.form.get('window')
    threshold = request.form.get('threshold')
    sequence_number = request.form.get('sequence_number')
    make_plots = True

    args = (0, chromosome, window, threshold, sequence_number, make_plots)  # 0 is a random index
    analysis_result, sequence_details, error = perform_analysis(args)
    if error:
        return error, 404

    analysis_summary_file = os.path.join(parent_dir,
                                         f'analysis_summary_chr_{chromosome}_window_{window}_threshold_{threshold}.csv')

    if os.path.exists(analysis_summary_file):
        df = pd.read_csv(analysis_summary_file)
        if int(sequence_number) not in df['Sequence_Number'].values:
            df = pd.concat([df, pd.DataFrame([sequence_details])], ignore_index=True)
            df.to_csv(analysis_summary_file, index=False)
    else:
        df = pd.DataFrame([sequence_details])
        df.to_csv(analysis_summary_file, index=False)

    return render_template('analyze.html', **analysis_result)


@app.route('/analyze_all', methods=['POST'])
def analyze_all():
    # Capture filter settings from the request
    # Similar to the index route's POST handling
    chr_selected = request.form.get('chromosome')
    window_selected = request.form.get('window')
    threshold_selected = request.form.get('threshold')
    alpha_min = float(request.form.get('alpha_min', 0))
    alpha_max = float(request.form.get('alpha_max', 1))
    beta_min = float(request.form.get('beta_min', 0))
    beta_max = float(request.form.get('beta_max', 1))
    mu_min = float(request.form.get('mu_min', 0))
    mu_max = float(request.form.get('mu_max', 1))
    invalid_min = int(request.form.get('invalid_min', 0))
    invalid_max = int(request.form.get('invalid_max', 100))
    start_index = int(request.form.get('start_index', 0)) - 1
    end_index = int(request.form.get('end_index', -1)) - 1  # Default to -1 to handle "0" as a valid end index

    if start_index < 0 or end_index < start_index:
        return "Invalid start or end index", 400

    folder_name = f"sequence_detection_{window_selected}_{threshold_selected}"
    file_path = os.path.join(CHR_DATA_PATH, chr_selected, folder_name, 'all_sequences.csv')

    if os.path.exists(file_path):
        df = pd.read_csv(file_path)
        filtered_data = df[(df['Alpha'] >= alpha_min) & (df['Alpha'] <= alpha_max) &
                           (df['Beta'] >= beta_min) & (df['Beta'] <= beta_max) &
                           (df['Mu'] >= mu_min) & (df['Mu'] <= mu_max) &
                           (df['Invalid'] >= invalid_min) & (df['Invalid'] <= invalid_max)]

        filtered_data = filtered_data.reset_index()
        filtered_data = filtered_data.iloc[start_index:end_index + 1]
        make_plots = False
        args_list = [(index, chr_selected, window_selected, threshold_selected, row['Sequence_number'], make_plots) for
                     index, row
                     in filtered_data.iterrows()]
        # Parallel processing
        pool = Pool(processes=os.cpu_count())  # Or a different number of processes
        results = pool.map(perform_analysis, args_list)
        pool.close()  # No more tasks
        pool.join()  # Wait for completion

        # Process results

        all_sequence_details = []
        for analysis_result, sequence_details, error in results:
            if not error:
                all_sequence_details.append(sequence_details)

        analysis_summary_file = os.path.join(parent_dir,
                                             f'analysis_summary_chr_{chr_selected}_window_{window_selected}_threshold_{threshold_selected}.csv')

        # Check if the file exists and append or create as necessary
        if os.path.exists(analysis_summary_file):
            existing_df = pd.read_csv(analysis_summary_file)
            new_df = pd.DataFrame(all_sequence_details)
            final_df = pd.concat([existing_df, new_df], ignore_index=True)
        else:
            final_df = pd.DataFrame(all_sequence_details)

        final_df.to_csv(analysis_summary_file, index=False)

        # Optionally, aggregate analysis_results or render them in a new template
        return render_template('analyze_all_results.html')
    else:
        return "Filtered data file not found", 404


@app.route('/analyze_all_chr', methods=['POST'])
def analyze_all_chr():
    window_selected = request.form.get('window')
    threshold_selected = request.form.get('threshold')
    alpha_min = float(request.form.get('alpha_min', 0))
    alpha_max = float(request.form.get('alpha_max', 1))
    beta_min = float(request.form.get('beta_min', 0))
    beta_max = float(request.form.get('beta_max', 1))
    mu_min = float(request.form.get('mu_min', 0))
    mu_max = float(request.form.get('mu_max', 1))
    invalid_min = int(request.form.get('invalid_min', 0))
    invalid_max = int(request.form.get('invalid_max', 100))

    folder_name = f"sequence_detection_{window_selected}_{threshold_selected}"
    chromosomes = [_ for _ in os.listdir(CHR_DATA_PATH) if os.path.isdir(os.path.join(CHR_DATA_PATH, _))]

    for chr_selected in chromosomes:
        file_path = os.path.join(CHR_DATA_PATH, chr_selected, folder_name, 'all_sequences.csv')

        if os.path.exists(file_path):
            df = pd.read_csv(file_path)
            filtered_data = df[(df['Alpha'] >= alpha_min) & (df['Alpha'] <= alpha_max) &
                               (df['Beta'] >= beta_min) & (df['Beta'] <= beta_max) &
                               (df['Mu'] >= mu_min) & (df['Mu'] <= mu_max) &
                               (df['Invalid'] >= invalid_min) & (df['Invalid'] <= invalid_max)]

            filtered_data = filtered_data.reset_index()
            make_plots = False
            args_list = [(index, chr_selected, window_selected, threshold_selected, row['Sequence_number'], make_plots)
                         for index, row in filtered_data.iterrows()]
            # Parallel processing
            pool = Pool(processes=os.cpu_count())  # Or a different number of processes
            results = pool.map(perform_analysis, args_list)
            pool.close()  # No more tasks
            pool.join()  # Wait for completion

            # Process results

            all_sequence_details = []
            for analysis_result, sequence_details, error in results:
                if not error:
                    all_sequence_details.append(sequence_details)

            analysis_summary_file = os.path.join(parent_dir,
                                                 f'analysis_summary_chr_{chr_selected}_window_{window_selected}_threshold_{threshold_selected}.csv')

            # Check if the file exists and append or create as necessary
            if os.path.exists(analysis_summary_file):
                existing_df = pd.read_csv(analysis_summary_file)
                new_df = pd.DataFrame(all_sequence_details)
                final_df = pd.concat([existing_df, new_df], ignore_index=True)
            else:
                final_df = pd.DataFrame(all_sequence_details)

            final_df.to_csv(analysis_summary_file, index=False)

    return render_template('analyze_all_results.html')


@app.route('/plot/<chromosome>/<window>/<threshold>/<sequence_number>/<decode_type>')
def serve_plot(chromosome, window, threshold, sequence_number, decode_type):
    file_path = os.path.join(CHR_DATA_PATH, chromosome,
                             f"sequence_detection_{window}_{threshold}",
                             "all_sequences", f"sequence_{int(float(sequence_number))}",
                             f"{decode_type}_decode", "plot.html")
    print(file_path)
    if os.path.exists(file_path):
        return send_file(file_path)
    else:
        return "File not found", 404


if __name__ == '__main__':
    app.run(host='0.0.0.0', port=8100, debug=False)
