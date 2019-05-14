"""
Implements the Kalman filter for a linear Gaussian state space model.
"""

import numpy as np
from numpy import dot
from scipy.linalg import inv
from textwrap import dedent

class Kalman:
    """
    Implements the Kalman filter for the Gaussian state space model
    .. math::
        x_{t+1} = A x_t + C w_{t+1} \\
        y_t = G x_t + H v_t
    Here :math:`x_t` is the hidden state and :math:`y_t` is the measurement.
    The shocks :math:`w_t` and :math:`v_t` are iid standard normals. Below
    we use the notation
    .. math::
        Q := CC'
        R := HH'
    Parameters
    -----------
    A : array_like or scalar(float)
        Part of the state transition equation.  It should be `n x n`
    C : array_like or scalar(float)
        Part of the state transition equation.  It should be `n x m`
    G : array_like or scalar(float)
        Part of the observation equation.  It should be `k x n`
    H : array_like or scalar(float), optional(default=None)
        Part of the observation equation.  It should be `k x l`
    x_hat : scalar(float) or array_like(float), optional(default=None)
        An n x 1 array representing the mean x_hat of the
        prior/predictive density.  Set to zero if not supplied.
    Sigma : scalar(float) or array_like(float), optional(default=None)
        An n x n array representing the covariance matrix Sigma of
        the prior/predictive density.  Must be positive definite.
        Set to the identity if not supplied.
    Attributes
    ----------
    Sigma, x_hat : as above
    """

    def __init__(self, A, C, G, H=None, x_hat=None, Sigma=None):
        self.A = A
        self.C = C
        self.G = G
        self.Q = np.dot(self.C, self.C.T)
        self.m = self.C.shape[1]
        self.k, self.n = self.G.shape

        if H is None:
            self.H = np.zeros((self.k, self.n))
        else:
            self.H = np.atleast_2d(H)
        if Sigma is None:
            self.Sigma = np.identity(self.n)
        else:
            self.Sigma = np.atleast_2d(Sigma)
        if x_hat is None:
            self.x_hat = np.zeros((self.n, 1))
        else:
            self.x_hat = np.atleast_2d(x_hat)
            self.x_hat.shape = self.n, 1

        self.R = np.dot(self.H, self.H.T)


    def __repr__(self):
        return self.__str__()

    def prior_to_filtered(self, y):
        r"""
        Updates the moments (x_hat, Sigma) of the time t prior to the
        time t filtering distribution, using current measurement :math:`y_t`.
        The updates are according to
        .. math::
            \hat{x}^F = \hat{x} + \Sigma G' (G \Sigma G' + R)^{-1}
                (y - G \hat{x})
            \Sigma^F = \Sigma - \Sigma G' (G \Sigma G' + R)^{-1} G
                \Sigma
        Parameters
        ----------
        y : scalar or array_like(float)
            The current measurement
        """
        # === simplify notation === #
        G, H = self.G, self.H
        R = np.dot(H, H.T)

        # === and then update === #
        y = np.atleast_2d(y)
        y.shape = self.k, 1
        E = dot(self.Sigma, G.T)
        F = dot(dot(G, self.Sigma), G.T) + R
        M = dot(E, inv(F))
        print(E.shape, F.shape, M.shape)
        self.x_hat = self.x_hat + dot(M, (y - dot(G, self.x_hat)))
        self.Sigma = self.Sigma - dot(M, dot(G, self.Sigma))

    def filtered_to_forecast(self):
        """
        Updates the moments of the time t filtering distribution to the
        moments of the predictive distribution, which becomes the time
        t+1 prior
        """
        # === simplify notation === #
        A, C = self.A, self.C
        Q = dot(C, C.T)

        # === and then update === #
        self.x_hat = dot(A, self.x_hat)
        self.Sigma = dot(A, dot(self.Sigma, A.T)) + Q

    def update(self, y):
        """
        Updates x_hat and Sigma given k x 1 nparray y.  The full
        update, from one period to the next
        Parameters
        ----------
        y : np.ndarray
            A k x 1 ndarray y representing the current measurement
        """
        self.prior_to_filtered(y)
        self.filtered_to_forecast()

    def __str__(self):
        m = """\
        Kalman filter:
          - dimension of state space          : {n}
          - dimension of observation equation : {k}
        """
        return dedent(m.format(n=self.n, k=self.k))

    def log_likelihood(self, y):
        """
        Computes log-likelihood of period ``t``
        Parameters
        ----------
        y : np.ndarray
            A k x 1 ndarray y representing the current measurement
        """

        eta = y - np.dot(self.G, self.x_hat)  # forecast error
        P = np.dot(self.G, np.dot(self.Sigma, self.G.T)) + self.R # covariance matrix of forecast error
        logL = - (y.shape[0] * np.log(2*np.pi) + np.log(np.linalg.det(P)) + np.sum(dot(eta.T, dot(inv(P), eta))))/2
        return logL

    def compute_loglikelihood(self, y):
        """
        Computes log-likelihood of entire observations
        Parameters
        ----------
        y : np.ndarray
            n x T matrix of observed data.
            n is the number of observed variables in one period.
           Each column is a vector of observations at each period.
        """
        T = y.shape[1]
        logL = 0
        # forecast and update
        for t in range(1, T):
            logL = logL + self.log_likelihood(y[:, t])
            self.update(y[:, t])
        return logL

