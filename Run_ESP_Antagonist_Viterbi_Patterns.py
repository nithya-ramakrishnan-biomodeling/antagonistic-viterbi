from Epigenetic_Sequence_Predictor.ESP_Sim_Viterbi import find_transition_point
import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import rcParams

# Find Transition Table and Heatmaps.
# mu_list = [0.0, 0.2, 0.3, 0.5, 0.8]
mu_list = [0.0, 0.2, 0.5, 0.8]
rho_list = [0.5]
# start_node_list = [1, 2]
# end_node_list = [1, 2]
start_node_list = [1]
end_node_list = [1]

alpha_list = np.round(np.arange(0.05, 1, 0.05), 2).tolist()
beta_list = np.round(np.arange(0.05, 1, 0.05), 2).tolist()
# For Extra-fine Uncomment the below.
# alpha_list = np.round(np.arange(0.01, 1, 0.01), 2).tolist()
# beta_list = np.round(np.arange(0.01, 1, 0.01), 2).tolist()

rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']
rcParams['font.size'] = 18
# rcParams['font.weight'] = 'bold'
rcParams['xtick.major.width'] = 2
rcParams['ytick.major.width'] = 2

for rho in rho_list:
    for mu in mu_list:
        for start_node in start_node_list:
            for end_node in end_node_list:
                print(f'Running for mu = {mu}, rho = {rho}, {start_node}x{end_node}')

                # Create a 2D list to hold transition points
                table = []
                for beta in beta_list:
                    row = []
                    for alpha in alpha_list:
                        k = find_transition_point(alpha, beta, mu=mu, rho=rho, start_node=start_node,
                                                  end_node=end_node)
                        if k is not None:
                            if k > 0:
                                k -= 1
                        else:
                            k = -100
                        row.append(k)
                    table.append(row)

                # Print the table
                # print(f"Transition Table for mu = {mu}:")
                # print("Beta/Alpha\t" + "\t".join([f"0.{i}" for i in range(1, 10)]))
                # for i, row in enumerate(table, start=1):
                #     print(f"0.{i}\t\t" + "\t".join([str(val) if val is not None else '-1' for val in row]))

                # Output the table to a CSV file
                with open(
                        f"transition_table_mu{str(mu).replace('.', '')}_rho{str(rho).replace('.', '')}_{start_node}x{end_node}.csv",
                        'w', newline='') as csvfile:
                    writer = csv.writer(csvfile)
                    import numpy as np

                    writer.writerow(['Beta/Alpha'] + [val for val in alpha_list])  # Header row
                    for i, row in enumerate(table, start=0):
                        writer.writerow([beta_list[i]] + [val if val is not None else -1 for val in row])

                # print("Transition table saved to transition_table.csv")

                # Define colors for fixed values and gradient
                fixed_colors = ['red', 'blue', 'green', 'yellow', 'grey']  # Colors for fixed values - adjust as needed

                # Create a figure
                fig, ax = plt.subplots(figsize=(7.5, 5.5), dpi=300)  # Adjust the size as needed (width, height)

                # Find the maximum value in the table
                max_val = np.max(table)

                # Loop through the table and plot squares with appropriate colors and values
                for i in range(len(table)):
                    for j in range(len(table[i])):
                        val = table[i][j]

                        # Set color for fixed values
                        if val == -500:
                            color = fixed_colors[0]
                            text_color = 'black'  # Text color for visibility
                        elif val == -300:
                            color = fixed_colors[1]
                            text_color = 'black'  # Text color for visibility
                        elif val == -200:
                            color = fixed_colors[2]
                            text_color = 'black'  # Text color for visibility
                        elif val == -100:
                            color = fixed_colors[3]
                            text_color = 'black'  # Text color for visibility
                        else:
                            # Calculate ratio to the maximum value for gradient color
                            ratio = (max_val - val) / max_val
                            color = (ratio, ratio, ratio)  # Shades of grey based on ratio
                            text_color = 'black'  # Text color for visibility

                        # Plot squares
                        rect = plt.Rectangle((j, i), 1, 1, facecolor=color, edgecolor='black')
                        ax.add_patch(rect)

                        # Add value text for grey region
                        show_k_values = True
                        if show_k_values:
                            if val > 0:
                                if val == max_val:
                                    # plt.text(j + 0.5, i + 0.5, str(val), color='white', ha='center', va='center')
                                    plt.text(j + 0.5, i + 0.5, str(val), color='white', ha='center', va='center',
                                             fontsize=12)
                                else:
                                    # plt.text(j + 0.5, i + 0.5, str(val), color=text_color, ha='center', va='center')
                                    plt.text(j + 0.5, i + 0.5, str(val), color=text_color, ha='center', va='center',
                                             fontsize=12)

                            else:
                                # plt.text(j + 0.5, i + 0.5, '-', color=text_color, ha='center', va='center')  # Show '-' for non-positive values
                                pass

                # Set the ticks and labels for x-axis and y-axis with rotation
                x_ticks = np.arange(0.5, len(table[0]), 1)
                y_ticks = np.arange(0.5, len(table), 1)

                # Create labels and make every other one an empty string for alternating active/inactive labels
                x_labels = [val if i % 2 == 0 else '' for i, val in enumerate(alpha_list)]
                y_labels = [val if i % 2 == 0 else '' for i, val in enumerate(beta_list)]

                ax.set_xticks(x_ticks)
                ax.set_xticklabels([val for val in x_labels], rotation=90)

                ax.set_yticks(y_ticks)
                ax.set_yticklabels([val for val in y_labels], rotation=0)

                # set custom tics and values for extrafine plots where we dont want all the values of alpha and beta
                # to be shown. Comment the default tick creation while uncomenting the bellow.
                # Find indices that match the desired steps of 0.1
                #
                # x_indices = np.arange(0.5, len(table[0]), 10)
                # x_values = alpha_list[0::10]
                # y_indices =  np.arange(0.5, len(table[0]), 10)
                # y_values = beta_list[0::10]
                #
                # # Set tick positions and labels for both axes
                # ax.set_xticks([x for x in x_indices])
                # ax.set_xticklabels([val for val in x_values], rotation=45, fontsize=12)
                #
                # ax.set_yticks([y for y in y_indices])
                # ax.set_yticklabels([val for val in y_values], rotation=45, fontsize=12)

                #######################################################################################################

                # Set labels and title
                # ax.set_title(f"Transition Heatmap for $\mu$ = {mu}, $\\rho$ = {rho}, ({start_node}x{end_node})")
                ax.set_xlabel("$\\alpha$")
                ax.set_ylabel("$\\beta$")

                # Set x and y axis limits to align the squares properly
                ax.set_xlim(0, len(table[0]))
                ax.set_ylim(0, len(table))

                # Create a legend for the colors
                legend_labels = ['-500: Error', '-300: All 1s', '-200: All 0s', '-100: Mix of 0s and 1s',
                                 '    K : Threshold value']

                legend_patches = [mpatches.Patch(color=fixed_colors[i], label=label) for i, label in
                                  enumerate(legend_labels)]

                # Add legend to the plot
                region_color = ['blue', 'grey', 'green', 'yellow']
                region_label = ['A', 'B', 'C', 'D']
                # region_color = ['blue', 'green', 'yellow']
                # region_label = ['A', 'C', 'D']
                legend_patches = [mpatches.Patch(color=region_color[i], label=label) for i, label in
                                  enumerate(region_label)]

                plt.legend(handles=legend_patches, loc='center left', bbox_to_anchor=(1, 0.85))

                # Save in all requested formats
                for fmt in ['png', 'tiff', 'eps']:
                    plt.savefig(
                        f"transition_heatmap_mu{str(mu).replace('.', '')}_rho{str(rho).replace('.', '')}_{start_node}x{end_node}.{fmt}",
                        bbox_inches='tight', dpi=300, format=fmt)
                plt.close()
                # Show plot
                # plt.show()
