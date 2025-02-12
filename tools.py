
import numpy as np

from scipy.sparse import diags, csr_array
from numpy.random import default_rng
from scipy.stats import norm

import matplotlib.pyplot as plt


def get_noise(b0, n_tot=1, gam0=1e-6, f_apx=1):
    """
    Simulates Poisson-Gaussian noise in signals.

    Parameters:
    b0      : array-like
              Raw signal representing counts or probability distribution function (PDF)
    n_tot   : float, optional
              Total counts over the PDF. Default is 1, treating b0 as counts.
    gam0    : float, optional
              Level of Gaussian background noise. Default is 1e6.
    f_apx   : bool, optional
              Flag to use Gaussian approximation for Poisson noise. Default is True.

    Returns:
    b       : array-like
              Data corrupted with Poisson-Gaussian noise
    Lb      : scipy.sparse matrix
              Cholesky factorization of the inverse of the covariance matrix
    """
    
    n_b = len(b0)  # Length of the data vector
    
    # POISSON NOISE
    theta = 1 / n_tot
    sig_pois = np.sqrt(theta * b0)
    
    # ADDITIVE GAUSSIAN NOISE
    gamma = np.max(sig_pois) * gam0
    sig_gaus = gamma
    
    # Combined noise standard deviation
    sig = np.sqrt(sig_pois**2 + sig_gaus**2)
    
    # Cholesky factorization of the inverse covariance matrix
    Lb = diags(1.0 / np.squeeze(sig), offsets=0, shape=(n_b, n_b), format='csr')
    
    # GENERATE NOISY DATA
    rng = default_rng(seed=0)  # Reset random number generator for consistent noise
    if f_apx == 1:  # Gaussian approximation
        epsilon = sig * rng.standard_normal(size=np.shape(b0))  # Gaussian noise
        b = np.maximum(np.round((b0 + epsilon) * n_tot), 0) / n_tot  # Add noise and remove negative counts
    else:  # Poisson noise
        b = sig_gaus * rng.standard_normal(size=np.shape(b0)) + rng.poisson(b0 * n_tot) / n_tot
        b = np.maximum(b, 0)  # Remove negative Gaussian noise
    
    return b, Lb


def gen_data(A, d, mu, s, w=None, d_star=None, N=1e3):
    """
    Generates data based on a multimodal normal distribution and adds noise.

    Parameters:
    A : array
        The matrix used for reconstruction.
    d : array
        Particle mobility diameters.
    mu : array
        Mean diameters for the distribution.
    s : array
        Standard deviations for the distribution.
    w : array, optional
        Weighting for each mode. Defaults to equal weights.
    d_star : array, optional
        Setpoints for plotting. If None, no plot is generated.
    N : int, optional
        Number of particles to simulate noise. Default is 1000.

    Returns:
    b : array
        Data with noise.
    Lb : array
        Noise covariance matrix.
    x0 : array
        The true underlying distribution.
    """

    # Weighting for multiple modes
    if w is None:
        w = np.ones(len(mu))

    # Determine whether to plot or not
    f_plot = d_star is not None

    # If N is not provided, set a default
    N = int(N)

    if not type(mu) == list:
        mu = [mu]
        s = [s]
        w = [w]

    # Initialize the true distribution x0
    x0 = np.zeros_like(d)

    # Sum over multiple modes (unimodal or multimodal distributions)
    for ii in range(len(mu)):
        x0 += w[ii] * norm.pdf(np.log(d), np.log(mu[ii]), np.log(s[ii])) * (np.sqrt(2 * np.pi) * np.log(s[ii]))

    # Calculate the noiseless signal b0
    b0 = A @ x0

    # Add noise to the data
    b, Lb = get_noise(b0, N)

    # Plot the data and distribution if d_star is provided
    if f_plot:
        d_star2 = np.exp(np.log(d_star) - (np.log(d_star[1]) - np.log(d_star[0])) / 2)  # Shift setpoints by half
        d_star2 = np.concatenate((d_star2, [np.exp(np.log(d_star[-1]) + (np.log(d_star[1]) - np.log(d_star[0])) / 2)]))

        plt.figure()
        plt.errorbar(d_star, b, 2 / Lb.diagonal(), fmt='none', ecolor=[0.8, 0.8, 0.8])

        plt.plot(d_star, b0, '-b', label='Real distribution')
        plt.stairs(b, d_star2, color=[0.8, 0.8, 0.8], label='Noisy data')
        plt.plot(d_star, b, 'b.', markersize=3)

        plt.xscale('log')
        plt.title('Real distribution and data')

        plt.ylim(bottom=0)
        plt.legend()
        plt.show()

    return b, Lb, x0
