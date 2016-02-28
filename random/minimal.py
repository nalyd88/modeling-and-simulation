# This script defines and tests a function for generating a sequence of
# uniformly distributed random numbers using the Minimal Standard generator
# developed by Park and Miller.
#
# Ref. 1: Numerical Recipes in C 2nd ed.
#

import matplotlib.pyplot as plt
import numpy as np


def min_rand(seed=0):
    """Minimal random number generator of Park and Miller. Returns a uniform
    random deviate between 0.0 and 1.0. Set or reset idum to any integer value
    (except the unlikely value MASK) to initialize the sequence; idum must not
    be altered between calls for successive deviates in a sequence.
    Ref. Numerical Recipes in C 2nd ed."""

    # Define constants
    ia = 16807
    im = 2147483647
    am = (1.0/im)
    iq = 127773
    ir = 2836
    mask = 123459876

    # Only allow the generator to be seeded once
    if "seed" not in min_rand.__dict__:
        # XORing with MASK allows use of zero and other simple bit patterns for
        # seed.
        min_rand.seed = seed ^ mask

    # Compute idum=(IA*idum) % IM without over-flows by Schrage's method.
    k = min_rand.seed/iq
    min_rand.seed = ia*(min_rand.seed - k*iq) - ir*k
    if min_rand.seed  < 0:
        min_rand.seed += im
    ans = am*min_rand.seed  # Convert to a floating result.
    min_rand.seed ^= mask   # Unmask before return.

    return ans


if __name__ == "__main__":

    # Test the distribution
    n = 10000
    nbins = int((2*n)**(1./3))
    x = [min_rand(0) for i in range(n)]
    plt.hist(x, bins=nbins, normed=True)
    plt.grid(True)
    plt.xlabel("Generated Numbers")
    plt.ylabel("Normalized Frequency (PDF)")
    plt.show()

    # Perform Chi-Squared Test
    E = 1
    O, bin_edges = np.histogram(x, bins=nbins, normed=True)
    Am = sum(O*E/E)
    print(Am, nbins)
