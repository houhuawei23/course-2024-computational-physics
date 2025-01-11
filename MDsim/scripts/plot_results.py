import os
import numpy as np
import matplotlib.pyplot as plt


def plot_results_beta(thermo_out_file: str, plot_file: str = "energy.png"):
    # Load data
    thermo = np.loadtxt(thermo_out_file)
    timeStep = 5 / 1000  # ps
    sampleInterval = 100
    timeInterval = timeStep * sampleInterval
    numData = thermo.shape[0]
    time = np.arange(1, numData + 1) * timeInterval
    totalEnergy = thermo[:, 1] + thermo[:, 2]
    relativeEnergy = (totalEnergy - np.mean(totalEnergy)) / np.abs(np.mean(totalEnergy))

    # Set up the figure and subplots
    plt.figure(figsize=(14, 10))
    plt.suptitle("Energy Analysis Over Time", fontsize=18, y=1.02)

    # Custom colors and styles
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]
    linewidth = 2.5

    # Plot 1: Kinetic Energy
    plt.subplot(2, 2, 1)
    plt.plot(time, thermo[:, 1], "-", color=colors[0], linewidth=linewidth)
    plt.xlabel("Time (ps)", fontsize=14)
    plt.ylabel("Kinetic Energy (eV)", fontsize=14)
    plt.title("(a) Kinetic Energy", fontsize=14)
    plt.grid(True, linestyle="--", alpha=0.6)

    # Plot 2: Potential Energy
    plt.subplot(2, 2, 2)
    plt.plot(time, thermo[:, 2], "-", color=colors[1], linewidth=linewidth)
    plt.xlabel("Time (ps)", fontsize=14)
    plt.ylabel("Potential Energy (eV)", fontsize=14)
    plt.title("(b) Potential Energy", fontsize=14)
    plt.grid(True, linestyle="--", alpha=0.6)

    # Plot 3: Total Energy
    plt.subplot(2, 2, 3)
    plt.plot(time, totalEnergy, "-", color=colors[2], linewidth=linewidth)
    plt.xlabel("Time (ps)", fontsize=14)
    plt.ylabel("Total Energy (eV)", fontsize=14)
    plt.title("(c) Total Energy", fontsize=14)
    plt.grid(True, linestyle="--", alpha=0.6)

    # Plot 4: Relative Energy Fluctuations
    plt.subplot(2, 2, 4)
    plt.plot(time[1:], relativeEnergy[1:], "-", color=colors[3], linewidth=linewidth)
    plt.xlabel("Time (ps)", fontsize=14)
    plt.ylabel("Relative Fluctuations", fontsize=14)
    plt.title("(d) Relative Energy Fluctuations", fontsize=14)
    plt.ylim((-1e-3, 1e-3))
    plt.grid(True, linestyle="--", alpha=0.6)

    # Adjust layout and save the figure
    plt.tight_layout()
    plt.savefig(plot_file, dpi=300, bbox_inches="tight")
    plt.show()


def check_energy_fluctuation_beta():
    # Define the directory path
    diff_cutoff_dir_path = os.path.join(os.getcwd(), "test/diff_cutoff")
    os.chdir(diff_cutoff_dir_path)

    # Load data from files
    energy_7A = np.sum(np.loadtxt("7A/thermo.out")[1:, 1:3], axis=1)
    energy_9A = np.sum(np.loadtxt("9A/thermo.out")[1:, 1:3], axis=1)
    energy_12A = np.sum(np.loadtxt("12A/thermo.out")[1:, 1:3], axis=1)
    energy_15A = np.sum(np.loadtxt("15A/thermo.out")[1:, 1:3], axis=1)

    # Calculate fluctuations
    fluctuation_7A = (energy_7A - np.mean(energy_7A)) / np.abs(np.mean(energy_7A))
    fluctuation_9A = (energy_9A - np.mean(energy_9A)) / np.abs(np.mean(energy_9A))
    fluctuation_12A = (energy_12A - np.mean(energy_12A)) / np.abs(np.mean(energy_12A))
    fluctuation_15A = (energy_15A - np.mean(energy_15A)) / np.abs(np.mean(energy_12A))

    # Prepare data for plotting
    x = ["7 A", "9 A", "12 A", "15 A"]
    y = [
        np.std(fluctuation_7A),
        np.std(fluctuation_9A),
        np.std(fluctuation_12A),
        np.std(fluctuation_15A),
    ]

    # Plotting
    plt.figure(figsize=(10, 8))
    bars = plt.bar(
        x, y, color=["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"], edgecolor="black"
    )

    # Add value labels on top of each bar
    for bar in bars:
        height = bar.get_height()
        plt.text(
            bar.get_x() + bar.get_width() / 2,
            height + 0.000001,  # Slightly above the bar
            f"{height:.6f}",
            ha="center",
            va="bottom",
            fontsize=12,
        )

    # Customize plot
    plt.ylabel("Relative Fluctuations", fontsize=14)
    plt.xlabel("Cutoff Distance", fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.title(
        "Comparison of Relative Fluctuations of Total Energy at Different Cutoffs",
        fontsize=16,
        pad=20,
    )
    plt.grid(axis="y", linestyle="--", alpha=0.7)

    # Save and show the plot
    # plt.tight_layout()
    plt.savefig("energy_fluctuation.png", dpi=300)
    plt.show()


if __name__ == "__main__":
    # thermo_out_file = './test/thermo.out'
    # plot_file = './test/energy.png'
    # plot_results_beta(thermo_out_file, plot_file)
    check_energy_fluctuation_beta()
