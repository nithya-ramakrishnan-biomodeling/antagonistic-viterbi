import pandas as pd
import matplotlib.pyplot as plt
import plotly.graph_objs as go
import os


def plot_chromosome_data_matplotlib(df, chromosome, output_dir):
    """
    Generates a plot using Matplotlib and saves it as a PNG file.
    """
    plt.figure(figsize=(12, 6))
    plt.plot(df['Alpha'], label='Alpha')
    plt.plot(df['Beta'], label='Beta')
    plt.plot(df['Mu'], label='Mu')

    # Adding markers for valid regions (where Invalid == 0)
    valid_indices = df[df['Invalid'] == 0].index
    plt.scatter(valid_indices, df.loc[valid_indices, 'Mu'], color='yellow', marker='o', label='Valid (Invalid=0)')

    plt.title(f"Matplotlib Plot: Sequence Data for {chromosome}")
    plt.xlabel("Sequence Number")
    plt.ylabel("Values")
    plt.legend()

    plt.savefig(os.path.join(output_dir, f"{chromosome}_matplotlib_sequence_plot.png"))
    plt.close()


def plot_chromosome_data_plotly(df, chromosome, output_dir):
    """
    Generates a plot using Plotly and saves it as an HTML file.
    """
    trace_alpha = go.Scatter(x=df.index, y=df['Alpha'], mode='lines', name='Alpha')
    trace_beta = go.Scatter(x=df.index, y=df['Beta'], mode='lines', name='Beta')
    trace_mu = go.Scatter(x=df.index, y=df['Mu'], mode='lines', name='Mu')

    # Adding markers for valid regions (where Invalid == 0)
    valid_indices = df[df['Invalid'] == 0].index
    marker_trace = go.Scatter(x=valid_indices, y=df.loc[valid_indices, 'Mu'], mode='markers',
                              marker=dict(color='yellow'), name='Valid (Invalid=0)')

    layout = go.Layout(
        title=f"Plotly Plot: Sequence Data for {chromosome}",
        xaxis=dict(title='Sequence Number', rangeslider=dict(visible=True)),
        yaxis=dict(title='Values'),
        showlegend=True
    )

    fig = go.Figure(data=[trace_alpha, trace_beta, trace_mu, marker_trace], layout=layout)
    fig.write_html(os.path.join(output_dir, f"{chromosome}_plotly_sequence_plot.html"))


if __name__ == "__main__":
    chr_data_folder_path = "Chr_Data"
    chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
                   "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
                   "chrX", "chrY", "chrMT"]

    for chromosome in chromosomes:
        print(f"Processing Chromosome :: {chromosome}")
        chromosome_folder_path = os.path.join(chr_data_folder_path, chromosome)

        if os.path.isdir(chromosome_folder_path):

            for folder_name in os.listdir(chromosome_folder_path):
                if folder_name.startswith("sequence_detection_"):
                    print(f"Working in {folder_name}")
                    detection_dir = os.path.join(chromosome_folder_path, folder_name)

                    # Process the all_sequences.csv file
                    all_sequences_file = os.path.join(detection_dir, "all_sequences.csv")
                    if os.path.exists(all_sequences_file):
                        df = pd.read_csv(all_sequences_file)

                        # Plot using Matplotlib
                        plot_chromosome_data_matplotlib(df, chromosome, detection_dir)

                        # Plot using Plotly
                        plot_chromosome_data_plotly(df, chromosome, detection_dir)
                    else:
                        print(f"All sequences file not found: {all_sequences_file}")

                    # Process each selected_sequences_{n} folder
                    for subfolder_name in os.listdir(detection_dir):
                        if subfolder_name.startswith("selected_sequences_"):
                            print(f"Working in {subfolder_name}")
                            selected_sequences_dir = os.path.join(detection_dir, subfolder_name)
                            valid_sequences_file = os.path.join(selected_sequences_dir, "valid_sequences.csv")

                            if os.path.exists(valid_sequences_file):
                                df = pd.read_csv(valid_sequences_file)

                                # Plot using Matplotlib
                                plot_chromosome_data_matplotlib(df, chromosome, selected_sequences_dir)

                                # Plot using Plotly
                                plot_chromosome_data_plotly(df, chromosome, selected_sequences_dir)
                            else:
                                print(f"Valid sequences file not found: {valid_sequences_file}")

    print("Completed!")
