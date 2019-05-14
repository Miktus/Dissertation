"""
Calculates the loss function to be further optimized
"""

from Solution_function_PyTorch import Solution
from Kalman_PyTorch import Kalman


def loss_function(params, F_initial, observation_matrix, observables, initial_Kalman):
    """
    Calculates the loss function to be further optimized

    Parameters
    -----------
    params : array_like or scalar(float)
        Parameters to be optimized for
    observation_matrix : array_like
        Array of indices of observable variables
    observables : array_like
        Array of observable variables
    """

    beta = params**2 / (1 + params**2)

    F, Q, nx, ny, nz = Solution(betta=beta, F_initial=F_initial)
    kalman = Kalman(A=F, C=Q, G=observation_matrix, x_hat=initial_Kalman)
    log_like = kalman.compute_loglikelihood(observables)

    return(log_like)
