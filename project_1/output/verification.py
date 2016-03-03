#!/usr/bin/env python
# This script provides functionality for model verification. More specifically,
# it provides functions for the analysis of output information produced by the
# main simulation.
# Dylan A. Crocker
# CSE 6730 Project 1
# Due: March 4, 2016

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm


def analyze_pedestrian_queue(file_path):
    """Analyze the pedestrians created by the simulation."""

    destinations, speeds = np.loadtxt(file_path, delimiter=',', unpack=True)

    # Show Histograms of destination selections
    plt.hist(destinations, bins=len(np.unique(destinations)), normed=True)
    plt.xlabel("Selected Destinations (IDs)")
    plt.ylabel("Normalized Frequency (PDF)")
    plt.grid(True)

    # Intended pedestrian speed distribution
    mu = 1.340
    sigma = 0.265
    pdf_range = np.arange(0, 2.5, 0.001)

    # Plot theoretical distribution
    plt.figure()
    plt.plot(pdf_range, norm.pdf(pdf_range, mu, sigma))
    plt.grid(True)
    plt.xlabel("Pedestrian Speeds (m/s)")
    plt.ylabel("PDF")

    # Plot the histogram of speeds with the theoretical distribution overlaid
    plt.figure()
    plt.hist(speeds, bins=int((2*len(speeds))**(1./3)), normed=True)
    plt.grid(True)
    plt.xlabel("Pedestrian Speeds (m/s)")
    plt.ylabel("Normalized Frequency (PDF)")
    plt.plot(pdf_range, norm.pdf(pdf_range, mu, sigma), linewidth=2)
    plt.legend(["Theory", "Actual"])

    plt.show()  # Show the plots


if __name__ == "__main__":

    # Define the files
    pedestrian_info_file = "pedestrians.txt"
    simulation_log_file = "log.txt"

    # Analyze files
    analyze_pedestrian_queue(pedestrian_info_file)
