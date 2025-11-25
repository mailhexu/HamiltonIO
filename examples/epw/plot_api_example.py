"""Example of using the plot_epw_distance API."""

import os
import matplotlib.pyplot as plt
from HamiltonIO.epw import plot_epw_distance


def main():
    # Find test data
    base = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
    data_dir = os.path.join(base, "HamiltonIO", "epw", "test", "up")

    # Create a figure and axes
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

    # Plot Rk-based distances
    plot_epw_distance(
        ax=ax1,
        path=data_dir,
        imode=0,
        distance_type="Rk",
        ylim=(1e-4, 10),
        color="blue",
    )
    ax1.set_title("Rk-based (WF-WF) distances")

    # Plot Rg-based distances
    plot_epw_distance(
        ax=ax2,
        path=data_dir,
        imode=0,
        distance_type="Rg",
        ylim=(1e-4, 10),
        color="green",
    )
    ax2.set_title("Rg-based (WF-atom) distances")

    plt.tight_layout()
    plt.savefig("epw_distance_plot_api_example.pdf")
    print("Plot saved to epw_distance_plot_api_example.pdf")


if __name__ == "__main__":
    main()
