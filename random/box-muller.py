# This script defines and tests a function for generating a sequence of normally
# distributed random numbers using the Box-Muller method.
#
# Ref. 1: Numerical Recipes in C 2nd ed.
# Ref. 2: Statistics for Engineers and Scientists 2nd ed. Navidi
#

import numpy as np
import matplotlib.pyplot as plt
import random


def gaussian(mu=0, sigma=1):
    """Generate a random number from a Gaussian distribution.

    Adapted from ref 1.
    """

    if "iset" not in gaussian.__dict__:
        gaussian.iset = 0
    if "gset" not in gaussian.__dict__:
        gaussian.gset = 0

    if gaussian.iset == 0:
        v1, v2, rsq = 0.0, 0.0, 0.0
        while rsq >= 1.0 or rsq == 0.0:
            v1 = 2.0*random.random() - 1
            v2 = 2.0*random.random() - 1
            rsq = v1*v1 + v2*v2
        fac = np.sqrt(-2.0*np.log(rsq)/rsq)
        gaussian.iset = 1
        gaussian.gset = v1*fac
        xi = v2*fac
    else:
        gaussian.iset = 0
        xi = gaussian.gset

    return mu + xi*sigma  # See ref 2 pg. 242


def gaussian_distribution(mu=0, sigma=1, n=100):
    """Create a list of random numbers from a Gaussian distribution."""
    return [gaussian(mu, sigma) for i in range(n)]


if __name__ == "__main__":

    # Test the distribution
    mean = 5
    stddev = 2
    nsamps = 100000
    nbins = int((2*nsamps)**(1./3))

    # Create the histogram of generated values
    x = gaussian_distribution(mu=mean, sigma=stddev, n=nsamps)
    plt.hist(x, bins=nbins, normed=True)
    plt.grid(True)
    plt.show()
