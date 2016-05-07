# Plot the simulation output

import numpy as np
import matplotlib.pyplot as plt


def plot_avg_hr(file_name):
    data = np.loadtxt(file_name, delimiter=',')
    plt.figure()
    plt.plot(data[:, 0], data[:, 1])
    plt.grid(True)
    plt.xlabel("Time (min)")
    plt.ylabel("Average HR")
    plt.show()


if __name__ == "__main__":
    plot_avg_hr("output_hr.txt")
