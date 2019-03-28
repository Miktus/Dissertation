% VERSION 2.0, MARCH 1997, COPYRIGHT H. UHLIG.
% EXAMPL0.M:
% Solving the stochastic neoclassical growth model with the "toolkit"
% Copyright: H. Uhlig.  Feel free to copy, modify and use at your own risk.
% However, you are not allowed to sell this software or otherwise impinge
% on its free distribution.
%--------------------------------------------------------------------------
% This version of exampl0 demonstrates the use of Martin Schneider's tool
% to create the system matrices AA-NN automatically from the equations
% of the model
% Usage: Copy the files callsysmat.m and makesysmat.m in the toolkit folder
%        and add 'callsysmat' as the first command in the file 'do_it.m'
%        Type 'help makesysmat' to learn more about the options
% Copyright: Martin Schneider, Oesterreichische Nationalbank, JUNE 2007
%--------------------------------------------------------------------------

clear;
clc;
close all

disp('DSGE');

%--------------------------------------------------------------------------
%-- Define the model
%--------------------------------------------------------------------------

%  VARIABLES   [M,3] cell array: Variable name (case sensitive) ~ Variable type ~ Description
%              Variable type: 'X' ... Endogenous state variable 
%                             'Y' ... Endogenous other (jump) variable  
%                             'Z' ... Exogenous state variable 
%                             'U' ... Innovation to exogenous state variable
%                             ''  ... Skip variable
VARIABLES = {...
            'e_a' 'Z' 'Aggregate productivity shock' 
            'e_I' 'Z' 'Adjustment cost shock' 
            'e_b' 'Z' 'Preference shock' 
            'e_L' 'Z' 'Labor supply shock' 
            'e_G' 'Z' 'Public spending shock' 
            'e_pi' 'Z' 'Inflation target' 
            'pi' 'X' 'Inflation rate' 
            'w' 'X' 'Real wages' 
            'K' 'X' 'Capital'
            'Q' 'X' 'Q-Tobin ratio' 
            'I' 'X' 'Investment' 
            'C' 'X' 'Consumption' 
            'R' 'X' 'Monetary policy interest rate'
            'r' 'Y' 'Capital rental rate'
            'L' 'Y' 'Labor'
            'Y' 'X' 'Output'
            'pi_f' 'X' 'Flexible inflation rate'
            'w_f' 'X' 'Flexible wages'
            'K_f' 'X' 'Flexible capital'
            'Q_f' 'X' 'Flexible Q-Tobin ratio'
            'I_f' 'X' 'Flexible investment' 
            'C_f' 'X' 'Flexible consumption' 
            'R_f' 'X' 'Flexible monetary policy interest rate'
            'r_f' 'Y' 'Flexible capital rental rate'
            'L_f' 'Y' 'Flexible labor'
            'Y_f' 'X' 'Flexible output'
            'ee_a' 'U' 'Aggregate productivity shock error' 
            'ee_I' 'U' 'Adjustment cost shock error' 
            'ee_b' 'U' 'Preference shock error' 
            'ee_L' 'U' 'Labor supply shock error' 
            'ee_G' 'U' 'Public spending shock error' 
            'ee_pi' 'U' 'Inflation target error' 
            'n_p' 'Z' 'Prices error' 
            'n_Q' 'Z' 'Tobin Q-ratio error' 
            'n_R' 'Z' 'Monetary policy interest rate error' 
            'n_w' 'Z' 'Wages error' 
            'one' 'X' 'Temporary variable'
            };

                
%    EQUATIONS   [N,3] cell array: Equation type ~ Equation name ~ equation
%                Equation type: 'D' ... Deterministic equation
%                               'E' ... Expectational equation 
%                               'S' ... Shock equation
%                               ''  ... Skip equation
EQUATIONS  = {...
        'E' '1: Consumption' ... 
            '0 = - C(t) + (h/(1 + h))*C(t-1) + (1/(1 + h))*C(t+1) - ((1 - h)/((1 + h)*sigma_c))*(R(t) - pi(t+1)) + ((1 - h)/((1 + h)*sigma_c))*(e_b(t)-e_b(t+1))'
        'E' '2:  Investment' ... 
            '0 = - I(t) + (1/(1 + betta))*I(t-1) + (betta/(1 + betta))*I(t+1) + (adj/(1 + betta))*Q(t) + (betta*e_I(t+1) - e_I(t))/(1 + betta)'
        'E' '3:  Q-equation' ... 
            '0 = - Q(t)  -R(t) + pi(t+1) + ((1 - tau)/(1 - tau + r_k))*Q(t+1) + (r_k/(1 - tau + r_k))*r(t+1) + n_Q(t)'
        'E' '4:  Capital' ... 
            '0 = - K(t) + (1 - tau)*K(t-1) + tau*I(t-1)'
        'E' '5:  Prices' ... 
            '0 = - pi(t) + (betta/(1 + betta*gamma_p))*pi(t+1) + (gamma_p/(1 + betta*gamma_p))*pi(t-1) + (((1 - betta*xi_p)*(1 - xi_p))/((1 + betta*gamma_p)*xi_p))*(alphaa*r(t) + (1 - alphaa)*w(t) - e_a(t) + n_p(t))'
        'E' '6:  Wages' ... 
            '0 = (-1 - ((1/(1 + betta))*((1 - betta*xi_w)*(1 - xi_w))/((1 + (1/lambda_w)*((1 + lambda_w)*sigma_L))*xi_w)))*w(t) + (betta/(1 + betta))*w(t+1) + (1/(1 + betta))*w(t-1)+ (betta/(1 + betta))*pi(t+1) - ((1 + betta*gamma_w)/(1 + betta))*pi(t) + (gamma_w/( 1 + betta))*pi(t-1) - ((1/(1 + betta))*((1 - betta*xi_w)*(1 - xi_w))/((1 + (1/lambda_w)*((1 + lambda_w)*sigma_L))*xi_w))*(- sigma_L*L(t) - (sigma_c/(1 - h))*(C(t) - h*C(t-1)) - e_L(t) - n_w(t))'
        'D' '7:  Labor demand' ... 
            '0 = - L(t) + -w(t) + (1 + psi)*r(t) + K(t-1)'
        'D' '8:  Production' ... 
            '0 = - Y(t) + phi*e_a(t) + phi*alphaa*K(t-1) + phi*alphaa*psi*r(t) + phi*(1 - alphaa)*L(t)'
        'E' '9:     Goods market' ... 
            '0 = - Y(t) + (1 - tau*k - g)*C(t) + tau*k*I(t) + e_G(t)'
        'E' '10: Monetary policy' ... 
            '0 = - R(t) + rho*R(t-1) + ((1 - rho)*r_pi - r_dpi)*pi(t-1) + (1 - rho)*(1 - r_pi)*e_pi(t) + Y(t)*((1-rho)*r_Y+r_dY) - Y_f(t)*((1-rho)*r_Y+r_dY) + r_dpi*pi(t) - r_dY*Y(t-1) - r_dY*Y_f(t-1) + n_R(t)'
        'E' '11:  Flexible Consumption' ... 
            '0 = - C_f(t)     + (h/(1 + h))*C_f(t-1) + (1/(1 + h))*C_f(t+1) - ((1 - h)/((1 + h)*sigma_c))*(R_f(t) - pi_f(t+1)) + ((1 - h)/((1 + h)*sigma_c))*(e_b(t)-e_b(t+1))'
        'E' '12:  Flexible Investment' ... 
            '0 = - I_f(t)     + (1/(1 + betta))*I_f(t-1) + (betta/(1 + betta))*I_f(t+1) + (adj/(1 + betta))*Q_f(t) + (betta*e_I(t+1) - e_I(t))/(1 + betta)'
        'E' '13:  Flexible Q-equation' ... 
           '0 = - Q_f(t) - R_f(t) + pi_f(t+1) + ((1 - tau)/(1 - tau + r_k))*Q_f(t+1) + (r_k/(1 - tau + r_k))*r_f(t+1)'
        'E' '14:  Flexible Capital' ... 
            '0 = - K_f(t)     + (1 - tau)*K_f(t-1) + tau*I_f(t-1)'
        'E' '15:  Flexible Prices' ... 
            '0       = alphaa*r_f(t)+(1-alphaa)*w_f(t) - e_a(t)'
        'E' '16:  Flexible Wages' ... 
            '0 = - w_f(t)     + sigma_L*L_f(t) + (sigma_c/(1 - h))*(C_f(t) - h*C_f(t-1)) - e_L(t)'
        'D' '17:  Flexible Labor demand' ... 
            '0 = - L_f(t)     + -w_f(t) + (1 + psi)*r_f(t) + K_f(t-1)'
        'D' '18:  Flexible Production' ... 
            '0 = - Y_f(t)     + phi*e_a(t) + phi*alphaa*K_f(t-1) + phi*alphaa*psi*r_f(t) + phi*(1 - alphaa)*L_f(t)'
        'E' '19:  Flexible Goods market' ... 
            '0 = - Y_f(t)     + (1 - tau*k - g)*C_f(t) + tau*k*I_f(t) + e_G(t)'
        'S' '20:   Aggregate productivity shock' ... 
            'e_I(t+1)     = rho_I*e_I(t) + ee_I(t+1)'
        'S' '21:   Adjustment cost shock' ... 
            'e_b(t+1)     = rho_b*e_b(t) + ee_b(t+1)'
        'S' '22:   Preference shock' ... 
            'e_L(t+1)     = rho_L*e_L(t) + ee_L(t+1)'
        'S' '23:   Labor supply shock' ... 
            'e_G(t+1)     = rho_G*e_G(t) + ee_G(t+1)'
        'S' '24:   Public spending shock' ... 
            'e_a(t+1)     = rho_a*e_a(t) + ee_a(t+1)'
        'S' '25:   Inflation target' ... 
            'e_pi(t+1)    = rho_pi*e_pi(t)+ ee_pi(t+1)'
        'E' '26: Temporary first equation' ...
             'pi_f(t) = 0*one(t)'
        'E' '27: Temporary second equation' ...
             'one(t) = 0*one(t-1)'
        'S' '26: Prices error expectation' ... 
            'n_p(t) = 0'  
        'S' '27: Tobin Q-ratio error expectation' ... 
            'n_Q(t) = 0'  
        'S' '28: Monetary policy interest rate error expectation' ... 
            'n_R(t) = 0'  
        'S' '29: Wages error expectation' ... 
            'n_w(t) = 0'  
            };

%--------------------------------------------------------------------------   
%--- Setting parameters and calculating the steady state
%--------------------------------------------------------------------------   

    betta       = 0.99;
    tau         = 0.025;
    alphaa      = 0.3;
    psi         = 1/0.169;
    gamma_p     = 0.469;
    gamma_w     = 0.763;
    lambda_w    = 0.5;
    xi_p        = 0.908;
    xi_w        = 0.737; 
    sigma_L     = 2.4;
    sigma_c     = 1.353;
    h           = 0.573;
    phi         = 1.408;
    adj         = 1/6.771;
    r_k         = (1/betta)-1+tau;
    k           = 8.8;
    g           = 0.18;
    r_dpi       = 0.14;
    r_Y         = 0.099;
    r_dY        = 0.159;
    rho         = 0.961; 
    r_pi        = 1.684; 
       
    
    rho_L       = 0.889;
    rho_a       = 0.823;
    rho_b       = 0.855;
    rho_G       = 0.949;
    rho_pi      = 0.924;
    rho_I       = 0.927;
  
    sd_R        = 0.081;
    sd_p        = 0.16;
    sd_w        = 0.289;
    sd_Q        = 0.604; 
 
%--------------------------------------------------------------------------   
%--- Setting the options and call the toolkit
%--------------------------------------------------------------------------   
SYSMAT_DISPDIAG     = 1;    % 1 ...  Display diagnostics for equations and variables and suggest equation/variable types
SYSMAT_OUTPUT       = 2;
DISPLAY_IMMEDIATELY = 1 ;   % Display warning messages immediately

do_it;
