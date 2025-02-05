"""
Function for solving DSGE models
Final version written by Michal Miktus, April 2019
"""

from torch import zeros, cat, eye
from Linear_Time_Iteration_PyTorch import Linear_Time_Iteration


def Solution(F_initial=zeros((31, 31)), betta=0.99, tau=0.025, alphaa=0.3, psi=1/0.169, gamma_p=0.469,
             gamma_w=0.763, lambda_w=0.5, xi_p=0.908, xi_w=0.737, sigma_L=2.4,
             sigma_c=1.353, h=0.573, phi=1.408, adj=1/6.771,
             k=8.8, g=0.18, r_dpi=0.14, r_Y=0.099, r_dY=0.159, rho=0.961, r_pi=1.684,
             rho_L=0.889, rho_a=0.823, rho_b=0.855, rho_G=0.949, rho_pi=0.924, rho_I=0.927,
             sd_R=0.081, sd_p=0.16, sd_w=0.289, sd_Q=0.604):

    r_k = (1 / betta) - 1 + tau

    # --------------------------------------------------------------------------
    # --- Matrix form in spirit of Uhlig
    # --------------------------------------------------------------------------

    # Structure of the model in Uhligs notation:
    #   0 = AA x(t) + BB x(t-1) + CC y(t) + DD z(t)
    #   0 = E_t [ FF x(t+1) + GG x(t) + HH x(t-1) + JJ y(t+1) + KK y(t) + LL z(t+1) + MM z(t)]
    #   z(t+1) = NN z(t) + epsilon(t+1),    with E_t [ epsilon(t+1) ] = 0
    #
    #   x(t)      : Endogenous state variables
    #   y(t)      : Endogenous other variables
    #   z(t)      : Exogenous state variables
    #   epsilon(t): Innovation to exogenous state variable

    AA = zeros((4, 17))
    BB = zeros((4, 17))
    CC = zeros((4, 4))
    DD = zeros((4, 10))
    FF = zeros((17, 17))
    GG = zeros((17, 17))
    HH = zeros((17, 17))
    JJ = zeros((17, 4))
    KK = zeros((17, 4))
    LL = zeros((17, 10))
    MM = zeros((17, 10))
    NN = zeros((10, 10))

    matrix_names = [AA, BB, CC, DD, FF, GG, HH, JJ, KK, LL, MM, NN]

    nx = AA.shape[1]
    ny = AA.shape[0]
    nz = DD.shape[1]

    # --- 1: Consumption

    GG[0, 5] = -1  # C(t)
    HH[0, 5] = (h/(1+h))  # C(t-1)
    FF[0, 5] = (1/(1+h))  # C(t+1)
    GG[0, 6] = -((1-h)/((1+h)*sigma_c))  # R(t)
    FF[0, 0] = ((1-h)/((1+h)*sigma_c))  # pi(t+1)
    MM[0, 2] = ((1-h)/((1+h)*sigma_c))  # e_b(t)
    LL[0, 2] = -((1-h)/((1+h)*sigma_c))   # e_b(t+1)

    # --- 2:  Investment

    GG[1, 4] = -1  # I(t)
    HH[1, 4] = (1/(1+betta))  # I(t-1)
    FF[1, 4] = (betta/(1+betta))  # I(t+1)
    GG[1, 3] = (adj/(1+betta))  # Q(t)
    LL[1, 1] = (1+betta)*betta  # e_I(t+1)
    MM[1, 1] = -(1+betta)  # e_I(t)

    # --- 3:  Q-equation

    GG[2, 3] = -1  # Q(t)
    GG[2, 6] = -1  # R(t)
    FF[2, 0] = 1  # pi(t+1)
    FF[2, 3] = ((1-tau)/(1-tau+r_k))  # Q(t+1)
    JJ[2, 0] = (r_k/(1-tau+r_k))  # r(t+1)
    MM[2, 7] = 1  # n_Q(t)

    # --- 4:  Capital

    GG[3, 2] = -1  # K(t)
    HH[3, 2] = (1-tau)  # K(t-1)
    HH[3, 4] = tau  # I(t-1)

    # --- 5:  Prices

    GG[4, 0] = -1  # pi(t)
    FF[4, 0] = (betta/(1+betta*gamma_p))  # pi(t+1)
    HH[4, 0] = (gamma_p/(1+betta*gamma_p))  # pi(t-1)
    KK[4, 0] = (((1-betta*xi_p)*(1-xi_p))/((1+betta*gamma_p)*xi_p))*alphaa  # r(t)
    GG[4, 1] = (((1-betta*xi_p)*(1-xi_p))/((1+betta*gamma_p)*xi_p))*(1-alphaa)  # w(t)
    MM[4, 0] = -(((1-betta*xi_p)*(1-xi_p))/((1+betta*gamma_p)*xi_p))  # e_a(t)
    MM[4, 6] = (((1-betta*xi_p)*(1-xi_p))/((1+betta*gamma_p)*xi_p))  # n_p(t)

    # --- 6:  Wages

    GG[5, 1] = (-1-((1/(1+betta))*((1-betta*xi_w)*(1-xi_w))/((1+(1/lambda_w)*((1+lambda_w)*sigma_L))*xi_w)))  # w(t)
    FF[5, 1] = (betta/(1+betta))  # w(t+1)
    HH[5, 1] = (1/(1+betta))  # w(t-1)
    FF[5, 0] = (betta/(1+betta))  # pi(t+1)
    GG[5, 0] = -((1+betta*gamma_w)/(1+betta))  # pi(t)
    HH[5, 0] = (gamma_w/(1+betta))   # pi(t-1)
    KK[
        5, 1] = ((1/(1+betta))*((1-betta*xi_w)*(1-xi_w))/((1+(1/lambda_w)*((1+lambda_w)*sigma_L))*xi_w))*sigma_L  # L(t)
    HH[
        5, 5] = -((1/(1+betta))*((1-betta*xi_w)*(1-xi_w))/((1+(1/lambda_w)*((1+lambda_w)*sigma_L))*xi_w))*(sigma_c/(1-h))*h  # C(t-1)
    GG[
        5, 5] = ((1/(1+betta))*((1-betta*xi_w)*(1-xi_w))/((1+(1/lambda_w)*((1+lambda_w)*sigma_L))*xi_w))*(sigma_c/(1-h))  # C(t)
    MM[5, 3] = ((1/(1+betta))*((1-betta*xi_w)*(1-xi_w))/((1+(1/lambda_w)*((1+lambda_w)*sigma_L))*xi_w))  # e_L(t)
    MM[5, 9] = ((1/(1+betta))*((1-betta*xi_w)*(1-xi_w))/((1+(1/lambda_w)*((1+lambda_w)*sigma_L))*xi_w))  # n_w(t)

    # --- 7:  Labor demand

    CC[0, 1] = -1  # L(t)
    AA[0, 1] = -1  # w(t)
    CC[0, 0] = (1+psi)  # r(t)
    BB[0, 2] = 1  # K(t-1)

    # --- 8:  Production

    AA[1, 7] = -1  # Y(t)
    DD[1, 0] = phi  # e_a(t)
    BB[1, 2] = phi*alphaa  # K(t-1)
    CC[1, 0] = phi*alphaa*psi  # r(t)
    CC[1, 1] = phi*(1-alphaa)  # L(t)

    # --- 9:     Goods market

    GG[6, 7] = -1  # Y(t)
    GG[6, 5] = (1-tau*k-g)  # C(t)
    GG[6, 4] = tau*k  # I(t)
    MM[6, 4] = 1  # e_G(t)

    # --- 10: Monetary policy

    GG[7, 6] = -1  # R(t)
    HH[7, 6] = rho  # R(t-1)
    HH[7, 0] = ((1-rho)*r_pi-r_dpi)  # pi(t-1)
    MM[7, 5] = (1-rho)*(1-r_pi)  # e_pi(t)
    GG[7, 7] = ((1-rho)*r_Y+r_dY)  # Y(t)
    GG[7, 15] = -((1-rho)*r_Y+r_dY)  # Y_f(t)
    GG[7, 0] = r_dpi  # pi(t)
    HH[7, 7] = -r_dY  # Y(t-1)
    HH[7, 15] = -r_dY  # Y_f(t-1)
    MM[7, 8] = 1  # n_R(t)

    # --- 11:  Flexible Consumption

    GG[8, 13] = -1  # C_f(t)
    HH[8, 13] = (h/(1+h))  # C_f(t-1)
    FF[8, 13] = (1/(1+h))  # C_f(t+1)
    GG[8, 14] = -((1-h)/((1+h)*sigma_c))  # R_f(t)
    FF[8, 8] = ((1-h)/((1+h)*sigma_c))  # pi_f(t+1)
    MM[8, 2] = ((1-h)/((1+h)*sigma_c))  # e_b(t)
    LL[8, 2] = -((1-h)/((1+h)*sigma_c))  # e_b(t+1)

    # --- 12:  Flexible Investment

    GG[9, 12] = -1  # I_f(t)
    HH[9, 12] = (1/(1+betta))  # I_f(t-1)
    FF[9, 12] = (betta/(1+betta))  # I_f(t+1)
    GG[9, 11] = (adj/(1+betta))  # Q_f(t)
    LL[9, 1] = (1+betta)*betta  # e_I(t+1)
    MM[9, 1] = -(1+betta)  # e_I(t)

    # --- 13:  Flexible Q-equation

    GG[10, 11] = -1  # Q_f(t)
    GG[10, 14] = -1  # R_f(t)
    FF[10, 8] = 1  # pi_f(t+1)
    FF[10, 11] = ((1-tau)/(1-tau+r_k))  # Q_f(t+1)
    JJ[10, 2] = (r_k/(1-tau+r_k))  # r_f(t+1)

    # --- 14:  Flexible Capital

    GG[11, 10] = -1  # K_f(t)
    HH[11, 10] = (1-tau)  # K_f(t-1)
    HH[11, 12] = tau  # I_f(t-1)

    # --- 15:  Flexible Prices

    KK[12, 2] = alphaa  # r_f(t)
    GG[12, 9] = (1-alphaa)  # w_f(t)
    MM[12, 0] = -1  # e_a(t)

    # --- 16:  Flexible Wages

    GG[13, 9] = -1  # w_f(t)
    KK[13, 3] = sigma_L  # L_f(t)
    HH[13, 13] = -(sigma_c/(1-h))*h  # C_f(t-1)
    GG[13, 13] = (sigma_c/(1-h))  # C_f(t)
    MM[13, 3] = -1  # e_L(t)

    # --- 17:  Flexible Labor demand

    CC[2, 3] = -1  # L_f(t)
    AA[2, 9] = -1  # w_f(t)
    CC[2, 2] = (1+psi)  # r_f(t)
    BB[2, 10] = 1  # K_f(t-1)

    # --- 18:  Flexible Production

    AA[3, 15] = -1  # Y_f(t)
    DD[3, 0] = phi  # e_a(t)
    BB[3, 10] = phi*alphaa  # K_f(t-1)
    CC[3, 2] = phi*alphaa*psi  # r_f(t)
    CC[3, 3] = phi*(1-alphaa)  # L_f(t)

    # --- 19:  Flexible Goods market

    GG[14, 15] = -1  # Y_f(t)
    GG[14, 13] = (1-tau*k-g)  # C_f(t)
    GG[14, 12] = tau*k  # I_f(t)
    MM[14, 4] = 1  # e_G(t)

    # --- 20:   Aggregate productivity shock

    NN[0, 1] = rho_I  # e_I(t)

    # --- 21:   Adjustment cost shock

    NN[1, 2] = rho_b  # e_b(t)

    # --- 22:   Preference shock

    NN[2, 3] = rho_L  # e_L(t)

    # --- 23:   Labor supply shock

    NN[3, 4] = rho_G  # e_G(t)

    # --- 24:   Public spending shock

    NN[4, 0] = rho_a  # e_a(t)

    # --- 25:   Inflation target

    NN[5, 5] = rho_pi  # e_pi(t)

    # --- 26: Temporary first equation

    GG[15, 8] = -1  # pi_f(t)
    GG[15, 16] = 0  # one(t)

    # --- 27: Temporary second equation

    GG[16, 16] = -1  # one(t)
    HH[16, 16] = 0  # one(t-1)

    # --- 26: Prices error expectation

    NN[6, 6] = 0  # n_p(t)

    # --- 27: Tobin Q-ratio error expectation

    NN[7, 7] = 0  # n_Q(t)

    # --- 28: Monetary policy interest rate error expectation

    NN[8, 8] = 0  # n_R(t)

    # --- 29: Wages error expectation

    NN[9, 9] = 0  # n_w(t)

    # --------------------------------------------------------------------------
    # --- Matrix form in spirit of Linear Time Iteration
    # --------------------------------------------------------------------------

    A = cat((cat((matrix_names[1], zeros((ny, ny)), zeros((ny, nz))), 1), cat((matrix_names[6], zeros((nx, ny)), zeros((nx, nz))), 1),
             cat((zeros((nz, nx)), zeros((nz, ny)), matrix_names[11]), 1)), 0)
    B = cat((cat((matrix_names[0], matrix_names[2], matrix_names[3]), 1), cat((matrix_names[5], matrix_names[8], matrix_names[10]), 1),
             cat((zeros((nz, nx)), zeros((nz, ny)), -eye((nz))), 1)), 0)
    C = cat((zeros((ny, nx + ny+nz)), cat((matrix_names[4], matrix_names[7], matrix_names[9]), 1), zeros((nz, nx+ny+nz))), 0)

    F_iter, Q_iter = Linear_Time_Iteration(A, B, C, F_initial, 1e-16, 1e-4)

    Q_iter = Q_iter[:, -10:]

    return F_iter, Q_iter, nx, ny, nz
