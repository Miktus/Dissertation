# Estimating DSGE by Maximum Likelihood in Python
# Author: Michal Miktus
# Date: 20.03.2019

# Import libraries

from numpy.testing import assert_allclose
import matplotlib.pyplot as plt
import seaborn as sn
from IPython.display import display
from statsmodels.tools.numdiff import approx_fprime, approx_fprime_cs
import statsmodels.api as sm
from pandas_datareader.data import DataReader
import pandas as pd
from scipy import optimize, signal
import numpy as np
from __future__ import division
%matplotlib inline

# Set printing options
np.set_printoptions(precision=3, suppress=True, linewidth=120)
pd.set_option('float_format', lambda x: '%.3g' % x, )

# Save the names of the equations, variables, and parameters
equation_names = [
    'Consumption', 'Investment', 'Q-equation', 'Capital', 'Prices', 'Wages', 'Labor demand', 'Production',
    'Goods market', 'Monetary policy', 'Flexible_Consumption', 'Flexible_Investment', 'Flexible_Q-equation',
    'Flexible_Capital', 'Flexible_Prices', 'Flexible_Wages', 'Flexible_Labor demand', 'Flexible_Production',
    'Flexible_Goods market', 'Flexible_Monetary policy'
    ]
len(equation_names)


variable_names = [
    'Aggregate productivity shock', 'Adjustment cost shock', 'Preference shock', 'Labor supply shock',
    'Public spending shock', 'Inflation target', 'Inflation rate', 'Real wages', 'Capital', 'Q-Tobin ratio',
    'Investment', 'Consumption', 'Monetary policy interest rate', 'Capital rental rate', 'Labor', 'Output', 'Flexible inflation rate', 'Flexible wages', 'Flexible capital', 'Flexible Q-Tobin ratio',
    'Flexible investment', 'Flexible consumption', 'Flexible monetary policy interest rate', 'Flexible capital rental rate', 'Flexible labor', 'Flexible output', 'Aggregate productivity shock error',
    'Adjustment cost shock error', 'Preference shock error', 'Labor supply shock error',
    'Public spending shock error', 'Inflation target error', 'Prices error',
    'Tobin Q-ratio error', 'Monetary policy interest rate error', 'Wages  error'
    ]

parameter_names = [
    'Habit persistance', 'Intertemporal elasticity of substitution', 'Discount factor',
    'Adjustment costs', 'Depreciation rate', 'Degree of price indexation',
    'Probability of price peg', 'Output elasticity to capital', 'Degree of wages indexation',
    'Probability of wages peg', 'Elasticity of labor', 'Labor Dixitz substitution parameter',
    'Inverse of the elasticity of the capital utilization cost function',
    'Proportion of fix costs of the production', 'Steady state proportions of capital on output',
    'Steady state proportions of public spending on output', 'Monetary policy weights',
    'Inflation gap weight', 'Output gap weight', 'Inflation gap difference weight',
    'Output gap difference weight', 'Revenue of the marginal utilization',
    'Aggregate productivity shock autoregressive part', 'Adjustment cost shock autoregressive part',
    'Preference shock autoregressive part', 'Labor supply shock autoregressive part',
    'Public spending shock autoregressive part', 'Inflation target autoregressive part',
    'Prices error standard deviation', 'Tobin Q-ratio error standard deviation',
    'Monetary policy interest rate error standard deviation', 'Wages  error standard deviation'
]

# Save symbolic forms for printing

variable_symbols = [r"e_a", r"e_I", r"e_b", r"e_L", r"e_G", r"e_pi", r"pi", r"w", r"K", r"Q",
                    r"I", r"C", r"R", r"r", r"L", r"Y", r"pi_f", r"w_f", r"K_f", r"Q_f", r"I_f", r"C_f",
                    r"R_f", r"r_f", r"L_f", r"Y_f", r"ee_a", r"ee_I", r"ee_b", r"ee_L", r"ee_G", r"ee_pi",
                    r"n_p", r"n_Q", r"n_R", r"n_w"]

contemporaneous_variable_symbols = [
    r"$%s_t$" % symbol for symbol in variable_symbols
]

lead_variable_symbols = [
    r"$%s_{t+1}$" % symbol for symbol in variable_symbols
]

parameter_symbols = [r"$h$", r"$\sigma_c$", r"$\beta$", r"$adj$", r"$\tau$", r"$\gamma_p$", r"$\xi_p$", r"$\alpha$",
                     r"$\gamma_w$", r"$\xi_w$", r"$\sigma_L$", r"$\lambda_w$", r"$\psi$", r"$\phi$", r"$k$", r"$g$",
                     r"$\rho$", r"$r_pi$", r"$r_Y$", r"$r_dpi$", r"$r_dY$", r"$r_k$", r"$\rho_a$", r"$\rho_I$",
                     r"$\rho_b$", r"$\rho_L$", r"$\rho_G$", r"$\rho_pi$", r"$sd_Q$", r"$sd_p$", r"$sd_w$", r"$sd_R$"]

# True parameters values

parameters = pd.DataFrame({
    'name': parameter_names,
    'value': [0.99, 0.025, 0.3, 1/0.169, .469, 0.763, 0.5, 0.908, 0.737, 2.4, 1.353, 0.573, 1.408, 1/6.771, (1/0.99)-1+0.025,
              8.8, 0.18, 0.14, 0.099, 0.159, 0.961, 1.684, 0.889, 0.823, 0.855, 0.949, 0.924, 0.927, 0.081, 0.16, 0.289, 0.604
              ]
})
parameters

# Main class


class DSGE(object):
    def __init__(self, params=None):
        # Model dimensions
        self.k_params = 31
        self.k_variables = 38

        # Initialize parameters
        if params is not None:
            self.update(params)

    def update(self, params):
        # Save deep parameters
        self.discount_rate = params[0]
        self.beta = params[0]
        self.tau = params[1]
        self.alpha = params[2]
        self.psi = params[3]
        self.gamma_p = params[4]
        self.gamma_w = params[5]
        self.lambda_w = params[6]
        self.xi_p = params[7]
        self.xi_w = params[8]
        self.sigma_L = params[9]
        self.sigma_c = params[10]
        self.h = params[11]
        self.phi = params[12]
        self.adj = params[13]
        self.r_k = params[14]
        self.k = params[15]
        self.g = params[16]
        self.r_dpi = params[17]
        self.r_Y = params[18]
        self.r_dY = params[19]
        self.rho = params[20]
        self.r_pi = params[21]
        self.rho_L = params[22]
        self.rho_a = params[23]
        self.rho_b = params[24]
        self.rho_G = params[25]
        self.rho_pi = params[26]
        self.rho_I = params[27]
        self.sd_R = params[28]
        self.sd_p = params[29]
        self.sd_w = params[30]
        self.sd_Q = params[31]


['e_a',
 'e_I',
 'e_b',
 'e_L',
 'e_G',
 'e_pi',
 'pi',
 'w',
 'K',
 'Q',
 'I',
 'C',
 'R',
 'r',
 'L',
 'Y',
 'pi_f',
 'w_f',
 'K_f',
 'Q_f',
 'I_f',
 'C_f',
 'R_f',
 'r_f',
 'L_f',
 'Y_f',
 'ee_a',
 'ee_I',
 'ee_b',
 'ee_L',
 'ee_G',
 'ee_pi',
 'n_p',
 'n_Q',
 'n_R',
 'n_w']
