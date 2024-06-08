import pandas as pd
import plotly.graph_objects as go
import os

# Define the list of chromosomes for which you want to process the data
chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
               "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
               "chrX", "chrY", "chrMT"]

# Initialize an empty DataFrame to hold concatenated data
concatenated_df = pd.DataFrame()

# Loop through each chromosome and concatenate the data
for n in chromosomes:
    file_path = f"Chr_Data/{n}/sequence_detection_5_05/all_sequences.csv"
    print(file_path)
    if os.path.exists(file_path):
        print('file_found')
        # Read the CSV file and append it to the concatenated DataFrame
        df = pd.read_csv(file_path)
        concatenated_df = pd.concat([concatenated_df, df], ignore_index=True)

# Ensure the DataFrame is not empty
if not concatenated_df.empty:
    # Plotting
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=concatenated_df['Sequence_number'], y=concatenated_df['Alpha'], mode='lines', name='Alpha'))
    fig.add_trace(go.Scatter(x=concatenated_df['Sequence_number'], y=concatenated_df['Beta'], mode='lines', name='Beta'))
    fig.add_trace(go.Scatter(x=concatenated_df['Sequence_number'], y=concatenated_df['Mu'], mode='lines', name='Mu'))

    # Update plot layout
    fig.update_layout(title='Alpha, Beta, Mu Curves', xaxis_title='Sequence Number', yaxis_title='Value')

    # Save the plot as an HTML file
    fig.write_html('all_sequences_plot.html')
    print("Plot saved as all_sequences_plot.html.")
else:
    print("No data available to plot.")
