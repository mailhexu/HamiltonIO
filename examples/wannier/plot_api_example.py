"""
Example of using the plot_wannier_distance API.

Run this script from the HamiltonIO root directory:
    python examples/wannier/plot_api_example.py

This creates a 2-panel plot showing Wannier hopping elements vs distance
for SrMnO3, demonstrating both scatter and binned visualizations.
"""

import os

import matplotlib.pyplot as plt

from HamiltonIO.wannier import plot_wannier_distance


def main():
    # Find example data directory
    base = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
    data_dir = os.path.join(base, "tests", "data", "SrMnO3_wannier")

    # Check if example data exists
    if not os.path.exists(os.path.join(data_dir, "SrMnO3_hr.dat")):
        print("Error: Example data not found.")
        print(f"Please ensure SrMnO3_hr.dat exists in {data_dir}")
        return

    # Create figure with single plot
    fig, ax = plt.subplots(figsize=(6, 5))

    # Plot using the API
    plot_wannier_distance(
        ax=ax,
        path=data_dir,
        prefix="SrMnO3",
        structure_file="SrMnO3.pwi",
        structure_format="espresso-in",
        ylim=(1e-4, 10),
        color="steelblue",
        s=3,
        alpha=0.6,
    )

    ax.set_title("Wannier Hopping vs Distance (SrMnO3)")

    plt.tight_layout()
    plt.savefig("wannier_distance_plot_api.pdf")
    print("Plot saved to wannier_distance_plot_api.pdf")


if __name__ == "__main__":
    main()
