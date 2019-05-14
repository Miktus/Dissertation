# Estimating DSGE by Maximum Likelihood in Python by PyTorch
# Author: Michal Miktus
# Date: 20.03.2019

# Useful to debug: with torch.autograd.set_detect_anomaly(True):

# Import libraries

from Solution_function_PyTorch import Solution
import torch
import time
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from torch import tensor, zeros, mm, randn, normal, cat, squeeze, unsqueeze, argmax, max
from torch.distributions import uniform
from Loss_function_PyTorch import loss_function

# Set printing options

sns.set()
np.set_printoptions(precision=3, suppress=True, linewidth=120)
torch.set_printoptions(precision=3, linewidth=120)
pd.set_option('float_format', lambda x: '%.3g' % x, )

# Declare path

path = "/Users/miktus/Documents/PSE/Dissertation/DSGE_estimation"

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
var_exo = variables.loc[variables['types'] == 'Z']

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

betta = tensor(0.99, requires_grad=True)
tau = tensor(0.025, requires_grad=True)
alphaa = tensor(0.3, requires_grad=True)
psi = tensor(1/0.169, requires_grad=True)
gamma_p = tensor(0.469, requires_grad=True)
gamma_w = tensor(0.763, requires_grad=True)
lambda_w = tensor(0.5, requires_grad=True)
xi_p = tensor(0.908, requires_grad=True)
xi_w = tensor(0.737, requires_grad=True)
sigma_L = tensor(2.4, requires_grad=True)
sigma_c = tensor(1.353, requires_grad=True)
h = tensor(0.573, requires_grad=True)
phi = tensor(1.408, requires_grad=True)
adj = tensor(1/6.771, requires_grad=True)
r_k = (1/betta)-1+tau
k = tensor(8.8, requires_grad=True)
g = tensor(0.18, requires_grad=True)
r_dpi = tensor(0.14, requires_grad=True)
r_Y = tensor(0.099, requires_grad=True)
r_dY = tensor(0.159, requires_grad=True)
rho = tensor(0.961, requires_grad=True)
r_pi = tensor(1.684, requires_grad=True)

rho_L = tensor(0.889, requires_grad=True)
rho_a = tensor(0.823, requires_grad=True)
rho_b = tensor(0.855, requires_grad=True)
rho_G = tensor(0.949, requires_grad=True)
rho_pi = tensor(0.924, requires_grad=True)
rho_I = tensor(0.927, requires_grad=True)

sd_R = tensor(0.081, requires_grad=True)
sd_p = tensor(0.16, requires_grad=True)
sd_w = tensor(0.289, requires_grad=True)
sd_Q = tensor(0.604, requires_grad=True)

# --------------------------------------------------------------------------
# --- Solution
# --------------------------------------------------------------------------

F, Q, nx, ny, nz = Solution(F_initial=zeros((31, 31)), betta=betta, tau=tau, alphaa=alphaa, psi=psi, gamma_p=gamma_p,
                            gamma_w=gamma_w, lambda_w=lambda_w, xi_p=xi_p, xi_w=xi_w, sigma_L=sigma_L,
                            sigma_c=sigma_c, h=h, phi=phi, adj=adj, k=k, g=g, r_dpi=r_dpi, r_Y=r_Y,
                            r_dY=r_dY, rho=rho, r_pi=r_pi, rho_L=rho_L, rho_a=rho_a, rho_b=rho_b, rho_G=rho_G,
                            rho_pi=rho_pi, rho_I=rho_I, sd_R=sd_R, sd_p=sd_p, sd_w=sd_w, sd_Q=sd_Q)

# --------------------------------------------------------------------------
# --- Monte Carlo design
# --------------------------------------------------------------------------

monte_carlo_iter = 1
par1_monte = zeros(monte_carlo_iter)
par2_monte = zeros(monte_carlo_iter)
par3_monte = zeros(monte_carlo_iter)
par4_monte = zeros(monte_carlo_iter)
par5_monte = zeros(monte_carlo_iter)
loss_monte = zeros(monte_carlo_iter)

for iter in range(monte_carlo_iter):

    start_time = time.time()

    print("Number of Monte Carlo experiment: " + str(iter))

    # --------------------------------------------------------------------------
    # --- Simulation
    # --------------------------------------------------------------------------

    T = 300  # Number of periods to simulate

    epsilon = tensor([normal(mean=0, std=sd_Q), normal(mean=0, std=sd_p), normal(mean=0, std=sd_w),
                      normal(mean=0, std=sd_R)])
    epsilon = cat((randn(6), epsilon))

    X_sim = zeros((nx+ny+nz, T))
    X_sim[:, 0] = squeeze(mm(Q, torch.t(unsqueeze(epsilon, 0))))

    for t in range(1, T):
        epsilon = tensor([normal(mean=0, std=sd_Q),
                          normal(mean=0, std=sd_p), normal(mean=0, std=sd_w),
                          normal(mean=0, std=sd_R)])
        epsilon = cat((randn(6), epsilon))
        X_sim[:, t] = squeeze(mm(F, torch.t(unsqueeze(X_sim[:, t-1].clone(), 0))) + mm(Q, torch.t(unsqueeze(epsilon, 0))))

    # Discard first 100 observations

    X_sim_discarded = X_sim[:, 100:]

    # Plot for consumption

    # plt.plot(X_sim[11,:].detach().numpy())
    # plt.show()

    # --------------------------------------------------------------------------
    # --- Estimation
    # --------------------------------------------------------------------------

    # Observables: investment, consumption, labor
    # Indexes: 10, 11, 14

    no = 3  # Number of observables

    observation_matrix = tensor([
        [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
        [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
        [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0.,
            0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.]], requires_grad=True)

    observables = X_sim_discarded[(10, 11, 14), :]

    # Plot likelihood function depending on one specific parameter value, eq. beta

    plot = False

    if plot:
        beta_to_graph = torch.arange(start=0., end=1.01, step=0.005, requires_grad=True)
        likelihood_to_graph = []
        likelihood_to_graph = [loss_function([beta, tensor(0.908), tensor(0.949), tensor(0.3), tensor(0.025)], F, observation_matrix, observables, X_sim_discarded[:, 0]) for beta in beta_to_graph]

        plt.plot(beta_to_graph.detach().numpy(), likelihood_to_graph, color='blue')
        plt.xlabel('Beta')
        plt.ylabel('Likelihood')
        plt.title('Likelihood depending on beta')
        # plt.show()
        plt.savefig(path + '/Graphs/Likelihood depending on beta.pdf')

    # Optimization

    optimization = True

    if optimization:

        number_of_tries = 1
        par1_final = zeros(number_of_tries)
        par2_final = zeros(number_of_tries)
        par3_final = zeros(number_of_tries)
        par4_final = zeros(number_of_tries)
        par5_final = zeros(number_of_tries)
        loss_final = zeros(number_of_tries)

        for j in range(number_of_tries):

            print("The try number: " + str(j))

            distribution = uniform.Uniform(torch.Tensor([3.]), torch.Tensor([10.]))
            par1 = distribution.sample(torch.Size([1, 1])).requires_grad_()
            par2 = distribution.sample(torch.Size([1, 1])).requires_grad_()
            par3 = distribution.sample(torch.Size([1, 1])).requires_grad_()
            distribution = uniform.Uniform(torch.Tensor([0.65]), torch.Tensor([1.]))
            par4 = distribution.sample(torch.Size([1, 1])).requires_grad_()
            distribution = uniform.Uniform(torch.Tensor([0.16]), torch.Tensor([0.22]))
            par5 = distribution.sample(torch.Size([1, 1])).requires_grad_()

            learning_rate = 1e-12
            n_iter = 100

            optimizer = torch.optim.SGD(params=[par1, par2, par3, par4, par5], lr=learning_rate)

            def closure():

                # Before the backward pass, use the optimizer object to zero all of the
                # gradients for the Tensors it will update (which are the learnable weights
                # of the model)
                optimizer.zero_grad()

                loss_value = loss_function([par1, par2, par3, par4, par5], F, observation_matrix, observables, X_sim_discarded[:, 0])

                # Backward pass: compute gradient of the loss with respect to model parameters

                loss_value.backward(retain_graph=True)

                return par1, par2, par3, par4, par5, loss_value

            # Calling the step function on an Optimizer makes an update to its parameters

            for i in range(n_iter):
                print("Optimization step number: " + str(i))
                par1_vector, par2_vector, par3_vector, par4_vector, par5_vector, loss_vector = optimizer.step(closure)

            par1_final[j] = par1_vector
            par2_final[j] = par2_vector
            par3_final[j] = par3_vector
            par4_final[j] = par4_vector
            par5_final[j] = par5_vector
            loss_final[j] = loss_vector

        loss_monte[iter] = max(loss_final)
        index = argmax(loss_final)
        par1_monte[iter] = par1_final[index]
        par2_monte[iter] = par2_final[index]
        par3_monte[iter] = par3_final[index]
        par4_monte[iter] = par4_final[index]
        par5_monte[iter] = par5_final[index]

        final = pd.DataFrame({
            'Beta': par1_monte.detach().numpy()**2/(1+par1_monte.detach().numpy()**2),
            'varsigma_p': par2_monte.detach().numpy()**2/(1+par2_monte.detach().numpy()**2),
            'rho_G': par3_monte.detach().numpy()**2/(1+par3_monte.detach().numpy()**2),
            'alphaa': par4_monte.detach().numpy()**2/(1+par4_monte.detach().numpy()**2),
            'tau': par5_monte.detach().numpy()**2/(1+par5_monte.detach().numpy()**2),
            'Likelihood': loss_monte.detach().numpy()
        })

        # Export to csv

        # final.to_csv(path + "/Results/MC_All_F_100_100.csv", index=False)

        print("--- %s seconds ---" % (time.time() - start_time))
