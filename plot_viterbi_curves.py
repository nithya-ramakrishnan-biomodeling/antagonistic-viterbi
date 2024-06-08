import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import rcParams


def straight_line_eqn(alpha, mu):
    return alpha / (2 - mu)


def parabola_eqn(alpha, mu):
    numerator = alpha ** 2
    denominator = 2 * (1 - alpha) * (1 - (mu / 2))

    # Handling division by zero using np.where
    # Replace division by zero cases with a large negative number (-10) to indicate error
    valid_denominator = np.where(denominator == 0, 1, denominator)  # Avoid division by zero
    beta_values = 1 - (numerator / valid_denominator)

    # Indicating where the division by zero would occur by setting beta to -10
    beta_values = np.where(denominator == 0, -10, beta_values)

    return beta_values


def draw_2d_plots(alpha, mu_list):
    plt.figure(figsize=(7.5, 5.5), dpi=300)
    colours = ['r', 'g', 'b', 'k']
    for index, mu in enumerate(mu_list):
        beta_parabola = parabola_eqn(alpha, mu)
        beta_straight = straight_line_eqn(alpha, mu)
        # plt.plot(alpha, beta_straight, f'{colours[index % 4]}', label=f'Straight line, mu = {mu}')
        # plt.plot(alpha, beta_parabola, f'{colours[index % 4]}', label=f'Parabola, mu = {mu}')
        plt.plot(alpha, beta_straight, f'{colours[index % 4]}')
        plt.plot(alpha, beta_parabola, f'{colours[index % 4]}', label=f'$\mu$ = {mu}')
    plt.ylim([0, 1])
    plt.xlabel(r'$\alpha$')
    plt.ylabel(r'$\beta$')
    # plt.title(r'Relationship between $\alpha$ and $\beta$ for different values of $\mu$')
    plt.legend(loc='lower left')

    # Add annotations
    annotations = {
        'C': (0.2, 0.6),
        'B': (0.7, 0.8),
        'D': (0.5, 0.1),
        'A': (0.8, 0.3)
    }
    for label, (x, y) in annotations.items():
        plt.annotate(label, xy=(x, y), xytext=(0, 0), textcoords="offset points",
                     ha='center', va='center',
                     bbox=dict(boxstyle="circle,pad=0.3", facecolor="white", edgecolor="black"))

    for fmt in ['png', 'tiff', 'eps']:
        plt.savefig(f'plot_viterbi_curves_2d.{fmt}', dpi=300, format=fmt)

    plt.close()


def filter_values(alpha, mu, beta):
    """
    Filters out points where beta is not within the 0 to 1 range.
    """
    mask = (beta >= 0) & (beta <= 1)
    return alpha[mask], mu[mask], beta[mask]


def draw_3d_plot(alpha, mu_values):
    fig = plt.figure(figsize=(14, 8))
    ax = fig.add_subplot(111, projection='3d')
    colormap = cm.plasma  # Choose a colormap here (e.g., viridis, plasma, inferno, magma)
    mu_normalized = (mu_values - mu_values.min()) / (mu_values.max() - mu_values.min())

    for i, mu in enumerate(mu_values):
        color = colormap(mu_normalized[i])  # Generate color based on normalized mu value
        beta_parabola = parabola_eqn(alpha, mu)
        beta_straight = straight_line_eqn(alpha, mu)
        alpha_parabola_filtered, mu_parabola_filtered, beta_parabola_filtered = filter_values(alpha,
                                                                                              np.full_like(alpha,
                                                                                                           mu),
                                                                                              beta_parabola)
        alpha_straight_filtered, mu_straight_filtered, beta_straight_filtered = filter_values(alpha,
                                                                                              np.full_like(alpha,
                                                                                                           mu),
                                                                                              beta_straight)

        # Plot with the same color for each mu
        ax.plot(alpha_straight_filtered, mu_straight_filtered, beta_straight_filtered, color=color)
        ax.plot(alpha_parabola_filtered, mu_parabola_filtered, beta_parabola_filtered, color=color)

    ax.set_xlabel('Alpha')
    ax.set_ylabel('Mu')
    ax.set_zlabel('Beta')
    ax.set_title('3D Visualization of Alpha, Beta, and Mu Relationships')
    ax.view_init(elev=20, azim=-35)
    plt.savefig('plot_viterbi_curves_3d.png')
    plt.show()
    plt.close()


rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']
rcParams['font.size'] = 18
# Set global font weight
# rcParams['font.weight'] = 'bold'

# Set global tick width
# rcParams['axes.linewidth'] = 2  # Set the thickness of the axes lines
rcParams['xtick.major.width'] = 2
rcParams['ytick.major.width'] = 2

alpha = np.round(np.arange(0, 1, 0.01), 2)
# mu_list = np.round(np.arange(0, 1, 0.2), 1)
mu_list = [0.0, 0.2, 0.5, 0.8]

draw_2d_plots(alpha, mu_list)

# mu_list = np.round(np.arange(0, 1, 0.01), 2)
# draw_3d_plot(alpha, mu_list)
