# Estimating DSGE by Maximum Likelihood in Python by PyTorch
# Author: Michal Miktus
# Date: 20.03.2019

# Useful to debug: with torch.autograd.set_detect_anomaly(True):

# Import libraries

import matplotlib.pyplot as plt
import seaborn as sn

import numpy as np
import pandas as pd
import torch
from pandas_datareader.data import DataReader
from torch import tensor, zeros, mm, randn, normal, cat, squeeze, unsqueeze
from Solution_function_PyTorch import Solution


# %matplotlib inline

# Set printing options

np.set_printoptions(precision=3, suppress=True, linewidth=120)
torch.set_printoptions(precision=3, linewidth=120)
pd.set_option('float_format', lambda x: '%.3g' % x, )

# --------------------------------------------------------------------------
# -- Define the model
# --------------------------------------------------------------------------

#  VARIABLES   [M,3] cell array: Variable name (case sensitive) ~ Variable type ~ Description
#              Variable type: 'X' ... Endogenous state variable
#                             'Y' ... Endogenous other (jump) variable
#                             'Z' ... Exogenous state variable
#                             'U' ... Innovation to exogenous state variable
#                             ''  ... Skip variable

variable_symbols = [r'$e_a$', r'$e_I$', r'$e_b$', r'$e_L$',
                    r'$e_G$', r'$e_pi$', r'$pi$', r'$w$', r'$K$',
                    r'$Q$', r'$I$', r'$C$', r'$R$', r'$r$', r'$L$',
                    r'$Y$', r'$pi_f$', r'$w_f$', r'$K_f$', r'$Q_f$',
                    r'$I_f$', r'$C_f$', r'$R_f$', r'$r_f$', r'$L_f$',
                    r'$Y_f$', r'$ee_a$', r'$ee_I$', r'$ee_b$', r'$ee_L$', r'$ee_G$', r'$ee_pi$', r'$n_p$', r'$n_Q$', r'$n_R$', r'$n_w$', r'$one$']

variable_types = ['Z',    'Z',    'Z',    'Z',    'Z',    'Z',
                  'X',    'X',    'X',    'X',    'X',    'X',
                  'X',    'Y',    'Y',    'X',    'X',    'X',
                  'X',    'X',    'X',    'X',    'X',    'Y',
                  'Y',    'X',    'U',    'U',    'U',    'U',
                  'U',    'U',    'Z',    'Z',    'Z',    'Z',
                  'X', ]

variable_names = ['Aggregate productivity shock', 'Adjustment cost shock',
                  'Preference shock', 'Labor supply shock',
                  'Public spending shock', 'Inflation target', 'Inflation rate', 'Real wages', 'Capital', 'Q-Tobin ratio', 'Investment',
                  'Consumption', 'Monetary policy interest rate', 'Capital rental rate', 'Labor', 'Output',
                  'Flexible inflation rate', 'Flexible wages', 'Flexible capital', 'Flexible Q-Tobin ratio',
                  'Flexible investment', 'Flexible consumption', 'Flexible monetary policy interest rate',
                  'Flexible capital rental rate', 'Flexible labor', 'Flexible output', 'Aggregate productivity shock error',
                  'Adjustment cost shock error', 'Preference shock error', 'Labor supply shock error', 'Public spending shock error',
                  'Inflation target error', 'Prices error', 'Tobin Q-ratio error', 'Monetary policy interest rate error', 'Wages error',
                  'Temporary variable']


variables = pd.DataFrame({
        'names': variable_names,
        'types': variable_types,
        'symbols': variable_symbols
    })

var_endo_states = variables.loc[variables['types'] == 'X']
var_endo_controls = variables.loc[variables['types'] == 'Y']
var_exo  = variables.loc[variables['types'] == 'Z']

#    EQUATIONS   [N,3] cell array: Equation type ~ Equation name ~ equation
#                Equation type: 'D' ... Deterministic equation
#                               'E' ... Expectational equation
#                               'S' ... Shock equation
#                               ''  ... Skip equation

equation_formulas = [
    '0 = - C(t) + (h/(1 + h))*C(t-1) + (1/(1 + h))*C(t+1) - ((1 - h)/((1 + h)*sigma_c))*(R(t) - pi(t+1)) + ((1 - h)/((1 + h)*sigma_c))*(e_b(t)-e_b(t+1))',
    '0 = - I(t) + (1/(1 + betta))*I(t-1) + (betta/(1 + betta))*I(t+1) + (adj/(1 + betta))*Q(t) + (betta*e_I(t+1) - e_I(t))/(1 + betta)',
    '0 = - Q(t)  -R(t) + pi(t+1) + ((1 - tau)/(1 - tau + r_k))*Q(t+1) + (r_k/(1 - tau + r_k))*r(t+1) + n_Q(t)',
    '0 = - K(t) + (1 - tau)*K(t-1) + tau*I(t-1)', '0 = - pi(t) + (betta/(1 + betta*gamma_p))*pi(t+1) + (gamma_p/(1 + betta*gamma_p))*pi(t-1) + (((1 - betta*xi_p)*(1 - xi_p))/((1 + betta*gamma_p)*xi_p))*(alphaa*r(t) + (1 - alphaa)*w(t) - e_a(t) + n_p(t))',
    '0 = (-1 - ((1/(1 + betta))*((1 - betta*xi_w)*(1 - xi_w))/((1 + (1/lambda_w)*((1 + lambda_w)*sigma_L))*xi_w)))*w(t) + (betta/(1 + betta))*w(t+1) + (1/(1 + betta))*w(t-1)+ (betta/(1 + betta))*pi(t+1) - ((1 + betta*gamma_w)/(1 + betta))*pi(t) + (gamma_w/( 1 + betta))*pi(t-1) - ((1/(1 + betta))*((1 - betta*xi_w)*(1 - xi_w))/((1 + (1/lambda_w)*((1 + lambda_w)*sigma_L))*xi_w))*(- sigma_L*L(t) - (sigma_c/(1 - h))*(C(t) - h*C(t-1)) - e_L(t) - n_w(t))', '0 = - L(t) + -w(t) + (1 + psi)*r(t) + K(t-1)',
    '0 = - Y(t) + phi*e_a(t) + phi*alphaa*K(t-1) + phi*alphaa*psi*r(t) + phi*(1 - alphaa)*L(t)',
    '0 = - Y(t) + (1 - tau*k - g)*C(t) + tau*k*I(t) + e_G(t)',
    '0 = - R(t) + rho*R(t-1) + ((1 - rho)*r_pi - r_dpi)*pi(t-1) + (1 - rho)*(1 - r_pi)*e_pi(t) + Y(t)*((1-rho)*r_Y+r_dY) - Y_f(t)*((1-rho)*r_Y+r_dY) + r_dpi*pi(t) - r_dY*Y(t-1) - r_dY*Y_f(t-1) + n_R(t)',
    '0 = - C_f(t) + (h/(1 + h))*C_f(t-1) + (1/(1 + h))*C_f(t+1) - ((1 - h)/((1 + h)*sigma_c))*(R_f(t) - pi_f(t+1)) + ((1 - h)/((1 + h)*sigma_c))*(e_b(t)-e_b(t+1))',
    '0 = - I_f(t) + (1/(1 + betta))*I_f(t-1) + (betta/(1 + betta))*I_f(t+1) + (adj/(1 + betta))*Q_f(t) + (betta*e_I(t+1) - e_I(t))/(1 + betta)',
    '0 = - Q_f(t) - R_f(t) + pi_f(t+1) + ((1 - tau)/(1 - tau + r_k))*Q_f(t+1) + (r_k/(1 - tau + r_k))*r_f(t+1)',
    '0 = - K_f(t) + (1 - tau)*K_f(t-1) + tau*I_f(t-1)',
    '0 = alphaa*r_f(t)+(1-alphaa)*w_f(t) - e_a(t)',
    '0 = - w_f(t) + sigma_L*L_f(t) + (sigma_c/(1 - h))*(C_f(t) - h*C_f(t-1)) - e_L(t)',
    '0 = - L_f(t) + -w_f(t) + (1 + psi)*r_f(t) + K_f(t-1)',
    '0 = - Y_f(t) + phi*e_a(t) + phi*alphaa*K_f(t-1) + phi*alphaa*psi*r_f(t) + phi*(1 - alphaa)*L_f(t)',
    '0 = - Y_f(t) + (1 - tau*k - g)*C_f(t) + tau*k*I_f(t) + e_G(t)',
    'e_I(t+1) = rho_I*e_I(t) + ee_I(t+1)',
    'e_b(t+1) = rho_b*e_b(t) + ee_b(t+1)',
    'e_L(t+1) = rho_L*e_L(t) + ee_L(t+1)',
    'e_G(t+1) = rho_G*e_G(t) + ee_G(t+1)',
    'e_a(t+1) = rho_a*e_a(t) + ee_a(t+1)',
    'e_pi(t+1) = rho_pi*e_pi(t)+ ee_pi(t+1)',
    'pi_f(t) = 0*one(t)',
    'one(t) = 0*one(t-1)',
    'n_p(t) = 0',
    'n_Q(t) = 0',
    'n_R(t) = 0',
    'n_w(t) = 0'
    ]

equation_type = ['E',    'E',    'E',    'E',    'E',    'E',
                 'D',    'D',    'E',    'E',    'E',    'E',
                 'E',    'E',    'E',    'E',    'D',    'D',
                 'E',    'S',    'S',    'S',    'S',    'S',
                 'S',    'E',    'E',    'S',    'S',    'S',
                 'S']

equation_names = ['Consumption', 'Investment', 'Q-equation', 'Capital', 'Prices', 'Wages', 'Labor demand', 'Production',
                  'Goods market', 'Monetary policy', 'Flexible Consumption', 'Flexible Investment', 'Flexible Q-equation',
                  'Flexible Capital', 'Flexible Prices', 'Flexible Wages', 'Flexible Labor demand', 'Flexible Production',
                  'Flexible Goods market', 'Aggregate productivity shock', 'Adjustment cost shock',
                  'Preference shock', 'Labor supply shock', 'Public spending shock', 'Inflation target', 'Temporary first equation', 'Temporary second equation', 'Prices error expectation', 'Tobin Q-ratio error expectation',
                  'Monetary policy interest rate error expectation', 'Wages error expectation']

equations = pd.DataFrame({
        'names': equation_names,
        'types': equation_type,
        'formulas': equation_formulas
    })

# --------------------------------------------------------------------------
# --- Setting parameters
# --------------------------------------------------------------------------

betta = tensor(0.99, requires_grad = True)
tau = tensor(0.025, requires_grad = True)
alphaa = tensor(0.3, requires_grad = True)
psi = tensor(1/0.169, requires_grad = True)
gamma_p = tensor(0.469, requires_grad = True)
gamma_w = tensor(0.763, requires_grad = True)
lambda_w = tensor(0.5, requires_grad = True)
xi_p = tensor(0.908, requires_grad = True)
xi_w = tensor(0.737, requires_grad = True)
sigma_L = tensor(2.4, requires_grad = True)
sigma_c = tensor(1.353, requires_grad = True)
h = tensor(0.573, requires_grad = True)
phi = tensor(1.408, requires_grad = True)
adj = tensor(1/6.771, requires_grad = True)
r_k = (1/betta)-1+tau
k = tensor(8.8, requires_grad = True)
g = tensor(0.18, requires_grad = True)
r_dpi = tensor(0.14, requires_grad = True)
r_Y = tensor(0.099, requires_grad = True)
r_dY = tensor(0.159, requires_grad = True)
rho = tensor(0.961, requires_grad = True)
r_pi = tensor(1.684, requires_grad = True)

rho_L = tensor(0.889, requires_grad = True)
rho_a = tensor(0.823, requires_grad = True)
rho_b = tensor(0.855, requires_grad = True)
rho_G = tensor(0.949, requires_grad = True)
rho_pi = tensor(0.924, requires_grad = True)
rho_I = tensor(0.927, requires_grad = True)

sd_R = tensor(0.081, requires_grad = True)
sd_p = tensor(0.16, requires_grad = True)
sd_w = tensor(0.289, requires_grad = True)
sd_Q = tensor(0.604, requires_grad = True)

# --------------------------------------------------------------------------
# --- Solution
# --------------------------------------------------------------------------

F, Q, nx, ny, nz = Solution(betta = betta, tau = tau, alphaa = alphaa, psi = psi, gamma_p = gamma_p,
             gamma_w = gamma_w, lambda_w = lambda_w, xi_p = xi_p, xi_w = xi_w, sigma_L = sigma_L,
             sigma_c = sigma_c, h = h, phi = phi, adj = adj, k = k, g = g, r_dpi = r_dpi, r_Y = r_Y,
             r_dY = r_dY, rho = rho, r_pi = r_pi, rho_L = rho_L, rho_a = rho_a, rho_b = rho_b, rho_G = rho_G,
             rho_pi = rho_pi, rho_I = rho_I, sd_R = sd_R, sd_p = sd_p, sd_w = sd_w, sd_Q = sd_Q)
#
# # --------------------------------------------------------------------------
# # --- Simulation
# # --------------------------------------------------------------------------
#
# T = 1000  # Number of periods to simulate
#
# epsilon = tensor([normal(mean = 0, std = sd_Q), normal(mean = 0, std = sd_p), normal(mean = 0, std = sd_w),
#                   normal(mean = 0, std = sd_R)])
# epsilon = cat((randn(6), epsilon))
# epsilon = cat((zeros((nx+ny)), epsilon))
#
# X_sim = zeros((nx+ny+nz, T))
# X_sim[:, 0] = squeeze(mm(Q, torch.t(unsqueeze(epsilon, 0))))
#
# for t in range(1, T):
#     epsilon = tensor([normal(mean=0, std=sd_Q),
#                         normal(mean=0, std=sd_p), normal(mean=0, std=sd_w),
#                         normal(mean=0, std=sd_R)])
#     epsilon = cat((randn(6), epsilon))
#     epsilon = cat((zeros((nx + ny)), epsilon))
#     X_sim[:, t] = squeeze(mm(F, torch.t(unsqueeze(X_sim[:, t-1].clone(), 0))) + mm(Q, torch.t(unsqueeze(epsilon, 0))))
#
# # Plot for consumption
#
# # plt.plot(X_sim[11,:].detach().numpy())
# # plt.show()
#
#
# # --------------------------------------------------------------------------
# # --- Estimation on real data (TO BE DONE)
# # --------------------------------------------------------------------------
# # Get some data
# # start='1984-01'
# # end = '2015-01'
# # labor = DataReader('HOANBS', 'fred', start=start, end=end)        # hours
# # consumption = DataReader('PCECC96', 'fred', start=start, end=end) # billions of dollars
# # investment = DataReader('GPDI', 'fred', start=start, end=end)     # billions of dollars
# # population = DataReader('CNP16OV', 'fred', start=start, end=end)  # thousands of persons
# # recessions = DataReader('USRECQ', 'fred', start=start, end=end)
# #
# # # Collect the raw values
# # raw = pd.concat((labor, consumption, investment, population.resample('QS').mean()), axis=1)
# # raw.columns = ['labor', 'consumption', 'investment', 'population']
# # raw['output'] = raw['consumption'] + raw['investment']
# #
# # # Make the data consistent with the model
# # y = np.log(raw.output * 10**(9-3) / raw.population)
# # n = np.log(raw.labor * (1e3 * 40) / raw.population)
# # c = np.log(raw.consumption * 10**(9-3) / raw.population)
# #
# # # Make the data stationary
# # y = y.diff()[1:]
# # n = n.diff()[1:]
# # c = c.diff()[1:]
# #
# # # Construct the final dataset
# # econ_observed = pd.concat((y, n, c), axis=1)
# # econ_observed.columns = ['output','labor','consumption']
# #
# # fig, ax = plt.subplots(figsize=(13,4))
# #
# # dates = econ_observed.index._mpl_repr()
# #
# # ax.plot(dates, econ_observed.output, label='Output')
# # ax.plot(dates, econ_observed.labor, label='Labor')
# # ax.plot(dates, econ_observed.consumption, label='Consumption')
# #
# # rec = recessions.resample('QS').last().loc[econ_observed.index[0]:].iloc[:, 0].values
# # ylim = ax.get_ylim()
# # ax.fill_between(dates, ylim[0]+1e-5, ylim[1]-1e-5, rec, facecolor='k', alpha=0.1)
# #
# # ax.xaxis.grid()
# # ax.legend(loc='lower left')
# # plt.show()
#
# # --------------------------------------------------------------------------
# # --- Estimation
# # --------------------------------------------------------------------------
#
# # Observables: labor, consumption, investment
# # Indexes: 14, 11, 10
#
# no = 3  # Number of observables
#
# observation_matrix =  tensor([
#     [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0.,
#         0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
#     [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0.,
#         0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
#     [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0.,
#         0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.]], requires_grad = True)
#
# # observation_matrix = zeros((no, nx+ny+nz), requires_grad=True)
# # observation_matrix[0,10] = observation_matrix[1,11] = observation_matrix[2,14] = 1
#
# observables = X_sim[(10,11,14), :]
#
# # Optimization
#
# from Loss_function_PyTorch import loss_function
#
# par = torch.rand(1, requires_grad=True)
# learning_rate = 1e-5
# n_iter = 50
#
# optimizer = torch.optim.SGD(params=[par], lr=learning_rate)
#
# # optimizer.zero_grad()
# # loss = loss_function(par, observation_matrix, observables,  X_sim[:, 0])
# # print(loss)
# # loss.backward()
# # optimizer.step()
#
# def closure():
#     # Before the backward pass, use the optimizer object to zero all of the
#     # gradients for the Tensors it will update (which are the learnable weights
#     # of the model)
#     optimizer.zero_grad()
#
#     # Without a constant term in the likelihood function:
#
#     loss_value = loss_function(par, observation_matrix, observables,  X_sim[:, 0])
#
#     # Backward pass: compute gradient of the loss with respect to model parameters
#
#     loss_value.backward(retain_graph=True)
#
#     return loss_value, par
#
# # Calling the step function on an Optimizer makes an update to its parameters
#
# loss_vector = torch.empty(n_iter)
# par_vector = torch.empty(n_iter)
#
# for i in range(n_iter):
#     print(i)
#     optimizer.step(closure)
#     loss_vector[i], par_vector[i] = closure()
#
# print(par_vector)
