import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import plotly.express as px
import plotly.graph_objects as go
from matplotlib import rcParams


def plot_mean_comparison(dataframe, summary_path, chromosome):
    # Calculate mean alpha and beta for each Mu
    mean_alpha_beta = dataframe.groupby('Mu')[['Alpha', 'Beta']].mean().reset_index()
    print(mean_alpha_beta)
    # Create a DataFrame to store the mean values for each condition
    mean_values = pd.DataFrame()

    # Define the conditions
    conditions = ['Avg_Bit_Error_No_Anta', 'Avg_Bit_Error_K_Fill_Antagonism', 'Avg_Bit_Error_K_Patch_Antagonism']

    # Calculate the mean for each condition and Mu
    for condition in conditions:
        mean_values[condition] = dataframe.groupby('Mu')[condition].mean()

    # Reset index to get Mu as a column
    mean_values.reset_index(inplace=True)

    # Construct plot save directory and path
    plot_dir = os.path.join(summary_path, 'plots', 'mean_comparison')
    os.makedirs(plot_dir, exist_ok=True)

    # Plot
    plt.figure(figsize=(7.5, 5.5), dpi=300)
    for condition in conditions:
        plt.plot(mean_values['Mu'].values, mean_values[condition].values, marker='o', linestyle='-', label=condition)

    plt.title(f'Mean Comparison Across Mu for {chromosome}')
    # plt.xlabel('Mu Values')
    plt.ylabel('Average Bit Error')
    plt.legend()

    # Annotate mean alpha and beta under each point
    for _, row in mean_alpha_beta.iterrows():
        mu = row['Mu']
        mean_alpha = row['Alpha']
        mean_beta = row['Beta']
        text = f'$\mu$={mu}\n$\\alpha$={mean_alpha:.4f}\n$\\beta$={mean_beta:.4f}'
        # Adjust the y position to place the text below the plot
        y_position = plt.gca().get_ylim()[0] - (plt.gca().get_ylim()[1] - plt.gca().get_ylim()[0]) * 0.07
        plt.text(mu, y_position, text, ha='center', va='top', fontsize=12)  # , rotation=45)

    plt.legend()
    plt.subplots_adjust(bottom=0.2)  # Adjust the bottom to make room for text annotations

    # Save plot
    for fmt in ['png', 'tiff', 'eps']:
        plot_save_path = os.path.join(plot_dir, f'{chromosome}_mean_comparison_line_plot.{fmt}')
        plt.savefig(plot_save_path, dpi=300, format=fmt)
    plt.close()

    print(f'Mean comparison line plot saved: {plot_save_path}')


def plot_violin_no_vs_antagonism(dataframe, summary_path, chromosome):
    # Calculate mean alpha and beta for each Mu
    mean_alpha_beta = dataframe.groupby('Mu')[['Alpha', 'Beta']].mean().reset_index()
    # Melt the DataFrame to have a single column for Avg Bit Error values and a new column for conditions
    dataframe_melted = dataframe.melt(id_vars=['Mu', 'Alpha', 'Beta'],
                                      value_vars=['Avg_Bit_Error_No_Anta', 'Avg_Bit_Error_K_Fill_Antagonism'],
                                      # ,'Avg_Bit_Error_K_Patch_Antagonism'],
                                      var_name='Condition', value_name='AvgBitError')

    # Map the long-form condition names to shorter, more manageable ones for plotting
    dataframe_melted['Condition'] = dataframe_melted['Condition'].map({
        'Avg_Bit_Error_No_Anta': 'No Antagonism',
        'Avg_Bit_Error_K_Fill_Antagonism': 'Antagonism',
        # 'Avg_Bit_Error_K_Patch_Antagonism': 'K_Patch Antagonism'
    })

    # Construct plot save directory and path
    plot_dir = os.path.join(summary_path, 'plots', 'no_vs_antagonism')
    os.makedirs(plot_dir, exist_ok=True)

    # Plot
    plt.figure(figsize=(7.5, 5.5), dpi=300)
    ax = sns.violinplot(x='Mu', y='AvgBitError', hue='Condition', data=dataframe_melted, split=False, palette='muted')
    # plt.title(f'Violin Plots of Average Bit Error by Mu Values - {chromosome}')
    plt.xlabel(r'$\mu$')
    plt.ylabel('Average BER')

    # Modify or remove the legend title
    legend = ax.legend()
    legend.set_title('')  # Set to an empty string to remove the title

    # # Annotate mean alpha, beta, and Mu under each violin using mean_alpha_beta directly
    # for _, row in mean_alpha_beta.iterrows():
    #     mu = row['Mu']
    #     mean_alpha = row['Alpha']
    #     mean_beta = row['Beta']
    #     text = f'$\mu$={mu}\n$\\alpha$={mean_alpha:.4f}\n$\\beta$={mean_beta:.4f}'
    #     # Find the position for the current Mu value
    #     pos = list(mean_alpha_beta['Mu']).index(mu)
    #     ax.text(pos, ax.get_ylim()[0] - (ax.get_ylim()[1] - ax.get_ylim()[0]) * 0.08, text, ha='center', va='top',
    #             fontsize=12)
    #
    # # Adjust layout to make room for the annotations
    # plt.subplots_adjust(bottom=0.2)

    # Save plot
    for fmt in ['png', 'tiff', 'eps']:
        plot_save_path = os.path.join(plot_dir, f'{chromosome}_no_vs_antagonism_violin_plot.{fmt}')
        plt.savefig(plot_save_path, dpi=300, format=fmt)
    plt.close()

    print(f'No Antagonism vs. Antagonism violin plot saved: {plot_save_path}')


def plot_violin_compare(dataframe, summary_path, chromosome):
    # Melt the DataFrame to have a single column for Delta Bit Error values and a new column for antagonism types
    dataframe_melted = dataframe.melt(id_vars=['Mu', 'Alpha', 'Beta'],
                                      value_vars=['Delta_Bit_Error_K_Fill_Antagonism',
                                                  'Delta_Bit_Error_K_Patch_Antagonism'],
                                      var_name='AntagonismType', value_name='DeltaBitError')

    # Map the long form antagonism type names to shorter, more manageable ones for plotting
    dataframe_melted['AntagonismType'] = dataframe_melted['AntagonismType'].map({
        'Delta_Bit_Error_K_Fill_Antagonism': 'K_Fill',
        'Delta_Bit_Error_K_Patch_Antagonism': 'K_Patch'
    })

    # Construct plot save directory and path
    plot_dir = os.path.join(summary_path, 'plots', 'compare_violin')
    os.makedirs(plot_dir, exist_ok=True)

    # Plot
    plt.figure(figsize=(7.5, 5.5), dpi=300)
    sns.violinplot(x='Mu', y='DeltaBitError', hue='AntagonismType', data=dataframe_melted, split=True)
    plt.title(f'Comparative Violin Plots of Delta Bit Error by Mu Values - {chromosome}')
    plt.xlabel(r'$\mu$')
    plt.ylabel('Delta Bit Error')

    # Save plot
    for fmt in ['png', 'tiff', 'eps']:
        plot_save_path = os.path.join(plot_dir, f'{chromosome}_compare_violin_plot.{fmt}')
        plt.savefig(plot_save_path, dpi=300, format=fmt)
    plt.close()

    print(f'Comparative violin plot saved: {plot_save_path}')


def plot_common_violin(dataframe, summary_path, chromosome, antagonism_type):
    # Construct plot save path
    plot_dir = os.path.join(summary_path, 'plots', antagonism_type, 'common_violin')
    os.makedirs(plot_dir, exist_ok=True)

    # Plot
    plt.figure(figsize=(7.5, 5.5), dpi=300)
    sns.violinplot(x='Mu', y=f'Delta_Bit_Error_{antagonism_type}', data=dataframe)
    # plt.title(f'Violin Plots of Delta Bit Error for {antagonism_type} by Mu Values - {chromosome}')
    plt.xlabel(r'$\mu$')
    plt.ylabel(f'$\Delta$ BER')

    # Save plot
    for fmt in ['png', 'tiff', 'eps']:
        plot_save_path = os.path.join(plot_dir, f'{chromosome}_common_violin_plot.{fmt}')
        plt.savefig(plot_save_path, dpi=300, format=fmt)
    plt.close()


def plot_group_violin(dataframe, summary_path, chromosome, antagonism_type):
    # Create directory for group plots
    group_plot_dir = os.path.join(summary_path, 'plots', antagonism_type, f'{chromosome}')
    os.makedirs(group_plot_dir, exist_ok=True)

    # Iterate and plot
    for index, (alpha, beta) in dataframe[['Alpha', 'Beta']].drop_duplicates().iterrows():
        subset_dataframe = dataframe[(dataframe['Alpha'] == alpha) & (dataframe['Beta'] == beta)]

        plt.figure(figsize=(7.5, 5.5), dpi=300)
        sns.violinplot(x='Mu', y=f'Delta_Bit_Error_{antagonism_type}', data=subset_dataframe)
        plt.title(f'$\\alpha$: {alpha}, $\\beta$: {beta}')
        plt.xlabel(r'$\mu$')
        plt.ylabel(f'Delta Bit Error {antagonism_type}')

        for fmt in ['png', 'tiff', 'eps']:
            plot_save_path = os.path.join(group_plot_dir, f'alpha_{alpha}_beta_{beta}_violin_plot.{fmt}')
            plt.savefig(plot_save_path, dpi=300, format=fmt)
        plt.close()


def plot_average_line(dataframe, summary_path, chromosome, antagonism_type):
    plot_dir = os.path.join(summary_path, 'plots', antagonism_type, 'average_line')
    os.makedirs(plot_dir, exist_ok=True)

    # Initialize a figure using Plotly for interactive plotting
    fig = px.line()

    # Iterate over each unique pair of Alpha and Beta
    for index, (alpha, beta) in dataframe[['Alpha', 'Beta']].drop_duplicates().iterrows():
        subset_dataframe = dataframe[(dataframe['Alpha'] == alpha) & (dataframe['Beta'] == beta)]

        # Group by Mu and calculate the mean
        avg_delta_error_by_mu = subset_dataframe.groupby('Mu')[
            f'Delta_Bit_Error_{antagonism_type}'].mean().reset_index()

        # Add trace to Plotly figure
        fig.add_scatter(x=avg_delta_error_by_mu['Mu'].values,
                        y=avg_delta_error_by_mu[f'Delta_Bit_Error_{antagonism_type}'].values,
                        mode='lines+markers', name=f'Alpha={alpha}, Beta={beta}')

    # Update layout
    fig.update_layout(title=f'Average Delta Bit Error by Mu for All Alpha-Beta Pairs - {chromosome}',
                      xaxis_title='Mu Values',
                      yaxis_title=f'Average Delta Bit Error {antagonism_type}',
                      legend_title='Alpha-Beta Pairs',
                      legend=dict(orientation="h"))

    # Determine save paths
    plot_html_save_path = os.path.join(plot_dir, f'{chromosome}_average_delta_error_line_plot.html')

    # Save HTML plot
    fig.write_html(plot_html_save_path)

    # For PNG, we'll use matplotlib to keep the original static plot saving functionality
    plt.figure(figsize=(7.5, 5.5), dpi=300)
    for index, (alpha, beta) in dataframe[['Alpha', 'Beta']].drop_duplicates().iterrows():
        subset_dataframe = dataframe[(dataframe['Alpha'] == alpha) & (dataframe['Beta'] == beta)]
        avg_delta_error_by_mu = subset_dataframe.groupby('Mu')[
            f'Delta_Bit_Error_{antagonism_type}'].mean().reset_index()
        plt.plot(avg_delta_error_by_mu['Mu'].values, avg_delta_error_by_mu[f'Delta_Bit_Error_{antagonism_type}'].values,
                 label=f'Alpha={alpha}, Beta={beta}', marker='o')
    plt.title(f'Average Delta Bit Error by Mu for All Alpha-Beta Pairs - {chromosome}')
    plt.xlabel(r'$\mu$')
    plt.ylabel(f'Average Delta Bit Error {antagonism_type}')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    plt.tight_layout()

    for fmt in ['png', 'tiff', 'eps']:
        plot_save_path = os.path.join(plot_dir, f'{chromosome}_average_delta_error_line_plot.{fmt}')
        plt.savefig(plot_save_path, dpi=300, format=fmt)
    plt.close()

    print(f'Line plot saved: {plot_save_path}')
    print(f'Interactive line plot saved: {plot_html_save_path}')


def plot_sequence_trends(dataframe, summary_path, chromosome):
    # Create the plot directory if it doesn't exist
    plot_dir = os.path.join(summary_path, 'plots', 'sequence_trend')
    os.makedirs(plot_dir, exist_ok=True)

    # Define the save path for the HTML file
    plot_save_path = os.path.join(plot_dir, f'{chromosome}_combined_sequence_trends.html')

    # Create a Plotly figure
    fig = go.Figure()

    # Add traces for Alpha, Beta, and Mu
    fig.add_trace(go.Scatter(x=dataframe.index, y=dataframe['Alpha'], mode='lines+markers', name='Alpha'))
    fig.add_trace(go.Scatter(x=dataframe.index, y=dataframe['Beta'], mode='lines+markers', name='Beta'))
    fig.add_trace(go.Scatter(x=dataframe.index, y=dataframe['Mu'], mode='lines+markers', name='Mu'))

    # Update layout
    fig.update_layout(
        title=f"Combined Sequence Trends for {chromosome}",
        xaxis_title="Sequence Index",
        yaxis_title="Value",
        legend_title="Parameter",
        height=600, width=1000,
        margin=dict(l=20, r=20, t=40, b=20)
    )

    # Optionally, add a range slider to the x-axis for zooming and panning through the sequence index
    fig.update_layout(xaxis=dict(
        rangeslider=dict(visible=True),
        type='linear'  # Change to 'category' if your index is categorical
    ))

    # Save the plot as an HTML file
    fig.write_html(plot_save_path)

    print(f'Combined sequence trends plot saved: {plot_save_path}')


if __name__ == '__main__':
    # summary_paths = ['Mouse_Data/Exp_Batch_Analysis/exp_100_1',
    #                  'Mouse_Data/Exp_Batch_Analysis/exp_100_2',
    #                  'Mouse_Data/Exp_Batch_Analysis/exp_100_3',
    #                  'Mouse_Data/Exp_Batch_Analysis/exp_100_4',
    #                  'Mouse_Data/Exp_Batch_Analysis/exp_100_inv_5',
    #                  'Mouse_Data/Exp_Batch_Analysis/exp_100_inv_10',
    #                  'Mouse_Data/Exp_Batch_Analysis/exp_100_inv_10_mu_07']
    # summary_paths = ['Mouse_Data/Exp_Batch_Analysis/exp_100_3', ]
    # summary_paths = ['Human_Data/Exp_Batch_Analysis/Exp_1']
    summary_paths = ['Mouse_Data/Exp_Batch_Analysis/exp_100_6' ]
    chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
                   "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19"]
    # "chrX"]  # , "chrY", "chrMT"] Y and MT have no data so commenting them out.

    window = '5'
    threshold = '05'

    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['Arial']
    rcParams['font.size'] = 12

    for summary_path in summary_paths:
        combined_df_list = []
        for chromosome in chromosomes:
            # Construct file path
            csv_file_path = os.path.join(summary_path,
                                         f'analysis_summary_chr_{chromosome}_window_{window}_threshold_{threshold}.csv')

            # Read the CSV file
            chromosome_data = pd.read_csv(csv_file_path)

            # Plot sequence trends with original, unrounded values
            plot_sequence_trends(chromosome_data, summary_path, chromosome)

            # Now we combine all the chromosomes to create a combined dataframe.
            csv_file_path = os.path.join(summary_path,
                                         f'analysis_summary_chr_{chromosome}_window_{window}_threshold_{threshold}.csv')
            if os.path.exists(csv_file_path):  # Check if the file exists to handle cases where it might not
                df = pd.read_csv(csv_file_path)
                df['Chromosome'] = chromosome  # Add chromosome column if you want to retain this info
                combined_df_list.append(df)

            # Round the values for Mu to 1 decimal places
            # chromosome_data['Mu'] = (chromosome_data['Mu'] * 20).round() / 20  # Rounding to 0.05s
            chromosome_data['Mu'] = chromosome_data['Mu'].round(1)
            chromosome_data['Alpha'] = chromosome_data['Alpha'].round(1)
            chromosome_data['Beta'] = chromosome_data['Beta'].round(1)

            for antagonism_type in ['K_Fill_Antagonism', 'K_Patch_Antagonism', 'NN']:
                plot_common_violin(chromosome_data, summary_path, chromosome, antagonism_type)
                plot_group_violin(chromosome_data, summary_path, chromosome, antagonism_type)
                plot_average_line(chromosome_data, summary_path, chromosome, antagonism_type)
            plot_violin_compare(chromosome_data, summary_path, chromosome)
            plot_violin_no_vs_antagonism(chromosome_data, summary_path, chromosome)
            plot_mean_comparison(chromosome_data, summary_path, chromosome)

        combined_df = pd.concat(combined_df_list, ignore_index=True)
        combined_csv_path = os.path.join(summary_path,
                                         f'analysis_summary_chr_ALL_window_{window}_threshold_{threshold}.csv')
        combined_df.to_csv(combined_csv_path, index=False)
        print(f'Combined chromosome data saved: {combined_csv_path}')

        # Construct file path
        csv_file_path = os.path.join(summary_path,
                                     f'analysis_summary_chr_ALL_window_{window}_threshold_{threshold}.csv')

        # Read the CSV file
        chromosome_data = pd.read_csv(csv_file_path)

        plot_sequence_trends(chromosome_data, summary_path, "ALL")

        # Round the values for Mu to 1 decimal places
        # chromosome_data['Mu'] = (chromosome_data['Mu'] * 20).round() / 20  # Rounding to 0.05s
        chromosome_data['Mu'] = chromosome_data['Mu'].round(1)
        chromosome_data['Alpha'] = chromosome_data['Alpha'].round(1)
        chromosome_data['Beta'] = chromosome_data['Beta'].round(1)

        # Now run the plots for the ALL Chromosome
        for antagonism_type in ['K_Fill_Antagonism', 'K_Patch_Antagonism', 'NN']:
            plot_common_violin(chromosome_data, summary_path, "ALL", antagonism_type)
            plot_group_violin(chromosome_data, summary_path, "ALL", antagonism_type)
            plot_average_line(chromosome_data, summary_path, "ALL", antagonism_type)
        plot_violin_compare(chromosome_data, summary_path, "ALL")
        plot_violin_no_vs_antagonism(chromosome_data, summary_path, "ALL")
        plot_mean_comparison(chromosome_data, summary_path, "ALL")
