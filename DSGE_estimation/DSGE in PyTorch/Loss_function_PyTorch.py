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

    beta = params[0]**2 / (1 + params[0]**2)
    xi_p = params[1]**2 / (1 + params[1]**2)
    rho_G = params[2]**2 / (1 + params[2]**2)
    alphaa = params[3]**2 / (1 + params[3]**2)
    tau = params[4]**2 / (1 + params[4]**2)

    F, Q, nx, ny, nz = Solution(betta=beta, xi_p=xi_p, tau=tau, alphaa=alphaa, rho_G=rho_G, F_initial=F_initial)
    kalman = Kalman(A=F, C=Q, G=observation_matrix, x_hat=initial_Kalman)
    log_like = kalman.compute_loglikelihood(observables)

    return(log_like)
