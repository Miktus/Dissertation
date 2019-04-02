"""
Class to estimate VAR by MLE
Final version written by Michal Miktus, April 2019
"""

import numpy as np
from scipy.optimize import minimize

np.random.seed(1)

class VAR:
    """
    **** VECTOR AUTOREGRESSION (VAR) MODELS ****
    ----------
    Parameters
    data : np.array
        Field to specify the time series data that will be used.
    lags : int
        Field to specify how many lag terms the model will have.
    integ : int (default : 0)
        Specifies how many time to difference the dependent variables.
    target : str (pd.DataFrame) or int (np.array) (default : None)
        By default, all columns will be selected as the dependent variables.
    """

    def __init__(self, data, lags, target=None, integ=0):

        # Latent Variables
        self.lags = lags
        self.integ = integ
        self.target = target
        self.model_name = "VAR(" + str(self.lags) + ")"

        # Format the dependant variables
        self.data = data

        # Difference data
        X = np.transpose(self.data)
        for order in np.arange(self.integ):
            X = np.asarray([np.diff(i) for i in X])
        self.data = X.T

        self.T = self.data.shape[0]
        self.ylen = self.data.shape[1]

        """
        Y : np.array
            Contains the length-adjusted time series (accounting for lags)
        """

        self.Y = (self.data[self.lags:, ]).T

    def _design(self):
        """ Creates a design matrix
        Z : np.array
        """

        Z = np.ones(((self.ylen * self.lags + 1), (self.T - self.lags)))

        row_count = 1
        for lag in np.arange(1, self.lags + 1):
            for reg in np.arange(self.ylen):
                Z[row_count, :] = self.data[:, reg][(self.lags - lag):-lag]
                row_count += 1

        return(Z)

    def OLS(self):
        """ Creates OLS coefficient matrix
        ----------
        Parameters:
        NULL
        ----------
        Returns
        The coefficient matrix B
        """

        Z = self._design()
        return np.dot(np.dot(self.Y, np.transpose(Z)), np.linalg.inv(np.dot(Z, np.transpose(Z))))

    def _neg_loglike(self, par):
        """ Calculate the MLE value, given the mean vector and variance matrix
        """

        Z = self._design()[1:]

        coef = np.reshape(par[0:self.lags * self.ylen**2],
                          (self.ylen, self.ylen * self.lags))
        coef_mean = par[self.lags * self.ylen
                        ** 2:self.lags * self.ylen**2 + self.ylen]
        coef_var = np.diag(par[self.lags * self.ylen**2 + self.ylen:])

        Y_0 = (self.Y.T - coef_mean).T
        Z_0 = (Z.T - np.tile(coef_mean, self.lags)).T

        logLik = -self.Y.shape[1] * self.ylen * np.log(2 * np.pi) * .5 - .5 * self.Y.shape[1] * np.log(np.abs(np.linalg.det(
            coef_var))) - .5 * np.trace(np.dot(np.dot((Y_0 - np.dot(coef, Z_0)).T, np.linalg.inv(coef_var)), Y_0 - np.dot(coef, Z_0)))

        return -logLik

    def MLE(self):
        """ Creates MLE coefficient matrix
        ----------
        Parameters:
        NULL
        ----------
        Returns
        The coefficient matrix MLE
        ----------
        It is based on the assumption of normality of errors
        """

        cons = []

        for i in np.arange(self.ylen):
            cons.append(
                dict({'type': 'ineq', 'fun': eval("lambda x: x[-" + str(i) + "]")}))

        # Make a list of initial parameter guesses

        initParams = np.repeat(
            1, self.lags * (self.ylen**2) + self.ylen + self.ylen)

        # Run the minimizer
        results = minimize(self._neg_loglike, initParams,
                           constraints=cons, method='COBYLA')

        # Print the results
        return(results.x)