%
% Status : main Dynare file
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

if isoctave || matlab_ver_less_than('8.6')
    clear all
else
    clearvars -global
    clear_persistent_variables(fileparts(which('dynare')), false)
end
tic0 = tic;
% Save empty dates and dseries objects in memory.
dates('initialize');
dseries('initialize');
% Define global variables.
global M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info estimation_info ys0_ ex0_
options_ = [];
M_.fname = 'SmetsWouters';
M_.dynare_version = '4.5.6';
oo_.dynare_version = '4.5.6';
options_.dynare_version = '4.5.6';
%
% Some global variables initialization
%
global_initialization;
diary off;
diary('SmetsWouters.log');
M_.exo_names = 'ee_a';
M_.exo_names_tex = 'ee\_a';
M_.exo_names_long = 'ee_a';
M_.exo_names = char(M_.exo_names, 'ee_I');
M_.exo_names_tex = char(M_.exo_names_tex, 'ee\_I');
M_.exo_names_long = char(M_.exo_names_long, 'ee_I');
M_.exo_names = char(M_.exo_names, 'ee_b');
M_.exo_names_tex = char(M_.exo_names_tex, 'ee\_b');
M_.exo_names_long = char(M_.exo_names_long, 'ee_b');
M_.exo_names = char(M_.exo_names, 'ee_L');
M_.exo_names_tex = char(M_.exo_names_tex, 'ee\_L');
M_.exo_names_long = char(M_.exo_names_long, 'ee_L');
M_.exo_names = char(M_.exo_names, 'ee_G');
M_.exo_names_tex = char(M_.exo_names_tex, 'ee\_G');
M_.exo_names_long = char(M_.exo_names_long, 'ee_G');
M_.exo_names = char(M_.exo_names, 'ee_pi');
M_.exo_names_tex = char(M_.exo_names_tex, 'ee\_pi');
M_.exo_names_long = char(M_.exo_names_long, 'ee_pi');
M_.exo_names = char(M_.exo_names, 'n_p');
M_.exo_names_tex = char(M_.exo_names_tex, 'n\_p');
M_.exo_names_long = char(M_.exo_names_long, 'n_p');
M_.exo_names = char(M_.exo_names, 'n_Q');
M_.exo_names_tex = char(M_.exo_names_tex, 'n\_Q');
M_.exo_names_long = char(M_.exo_names_long, 'n_Q');
M_.exo_names = char(M_.exo_names, 'n_R');
M_.exo_names_tex = char(M_.exo_names_tex, 'n\_R');
M_.exo_names_long = char(M_.exo_names_long, 'n_R');
M_.exo_names = char(M_.exo_names, 'n_w');
M_.exo_names_tex = char(M_.exo_names_tex, 'n\_w');
M_.exo_names_long = char(M_.exo_names_long, 'n_w');
M_.endo_names = 'e_a';
M_.endo_names_tex = 'e\_a';
M_.endo_names_long = 'e_a';
M_.endo_names = char(M_.endo_names, 'e_I');
M_.endo_names_tex = char(M_.endo_names_tex, 'e\_I');
M_.endo_names_long = char(M_.endo_names_long, 'e_I');
M_.endo_names = char(M_.endo_names, 'e_b');
M_.endo_names_tex = char(M_.endo_names_tex, 'e\_b');
M_.endo_names_long = char(M_.endo_names_long, 'e_b');
M_.endo_names = char(M_.endo_names, 'e_L');
M_.endo_names_tex = char(M_.endo_names_tex, 'e\_L');
M_.endo_names_long = char(M_.endo_names_long, 'e_L');
M_.endo_names = char(M_.endo_names, 'e_G');
M_.endo_names_tex = char(M_.endo_names_tex, 'e\_G');
M_.endo_names_long = char(M_.endo_names_long, 'e_G');
M_.endo_names = char(M_.endo_names, 'e_pi');
M_.endo_names_tex = char(M_.endo_names_tex, 'e\_pi');
M_.endo_names_long = char(M_.endo_names_long, 'e_pi');
M_.endo_names = char(M_.endo_names, 'one');
M_.endo_names_tex = char(M_.endo_names_tex, 'one');
M_.endo_names_long = char(M_.endo_names_long, 'one');
M_.endo_names = char(M_.endo_names, 'pi');
M_.endo_names_tex = char(M_.endo_names_tex, 'pi');
M_.endo_names_long = char(M_.endo_names_long, 'pi');
M_.endo_names = char(M_.endo_names, 'w');
M_.endo_names_tex = char(M_.endo_names_tex, 'w');
M_.endo_names_long = char(M_.endo_names_long, 'w');
M_.endo_names = char(M_.endo_names, 'K');
M_.endo_names_tex = char(M_.endo_names_tex, 'K');
M_.endo_names_long = char(M_.endo_names_long, 'K');
M_.endo_names = char(M_.endo_names, 'Q');
M_.endo_names_tex = char(M_.endo_names_tex, 'Q');
M_.endo_names_long = char(M_.endo_names_long, 'Q');
M_.endo_names = char(M_.endo_names, 'I');
M_.endo_names_tex = char(M_.endo_names_tex, 'I');
M_.endo_names_long = char(M_.endo_names_long, 'I');
M_.endo_names = char(M_.endo_names, 'C');
M_.endo_names_tex = char(M_.endo_names_tex, 'C');
M_.endo_names_long = char(M_.endo_names_long, 'C');
M_.endo_names = char(M_.endo_names, 'R');
M_.endo_names_tex = char(M_.endo_names_tex, 'R');
M_.endo_names_long = char(M_.endo_names_long, 'R');
M_.endo_names = char(M_.endo_names, 'r');
M_.endo_names_tex = char(M_.endo_names_tex, 'r');
M_.endo_names_long = char(M_.endo_names_long, 'r');
M_.endo_names = char(M_.endo_names, 'L');
M_.endo_names_tex = char(M_.endo_names_tex, 'L');
M_.endo_names_long = char(M_.endo_names_long, 'L');
M_.endo_names = char(M_.endo_names, 'Y');
M_.endo_names_tex = char(M_.endo_names_tex, 'Y');
M_.endo_names_long = char(M_.endo_names_long, 'Y');
M_.endo_names = char(M_.endo_names, 'pi_f');
M_.endo_names_tex = char(M_.endo_names_tex, 'pi\_f');
M_.endo_names_long = char(M_.endo_names_long, 'pi_f');
M_.endo_names = char(M_.endo_names, 'w_f');
M_.endo_names_tex = char(M_.endo_names_tex, 'w\_f');
M_.endo_names_long = char(M_.endo_names_long, 'w_f');
M_.endo_names = char(M_.endo_names, 'K_f');
M_.endo_names_tex = char(M_.endo_names_tex, 'K\_f');
M_.endo_names_long = char(M_.endo_names_long, 'K_f');
M_.endo_names = char(M_.endo_names, 'Q_f');
M_.endo_names_tex = char(M_.endo_names_tex, 'Q\_f');
M_.endo_names_long = char(M_.endo_names_long, 'Q_f');
M_.endo_names = char(M_.endo_names, 'I_f');
M_.endo_names_tex = char(M_.endo_names_tex, 'I\_f');
M_.endo_names_long = char(M_.endo_names_long, 'I_f');
M_.endo_names = char(M_.endo_names, 'C_f');
M_.endo_names_tex = char(M_.endo_names_tex, 'C\_f');
M_.endo_names_long = char(M_.endo_names_long, 'C_f');
M_.endo_names = char(M_.endo_names, 'R_f');
M_.endo_names_tex = char(M_.endo_names_tex, 'R\_f');
M_.endo_names_long = char(M_.endo_names_long, 'R_f');
M_.endo_names = char(M_.endo_names, 'r_f');
M_.endo_names_tex = char(M_.endo_names_tex, 'r\_f');
M_.endo_names_long = char(M_.endo_names_long, 'r_f');
M_.endo_names = char(M_.endo_names, 'L_f');
M_.endo_names_tex = char(M_.endo_names_tex, 'L\_f');
M_.endo_names_long = char(M_.endo_names_long, 'L_f');
M_.endo_names = char(M_.endo_names, 'Y_f');
M_.endo_names_tex = char(M_.endo_names_tex, 'Y\_f');
M_.endo_names_long = char(M_.endo_names_long, 'Y_f');
M_.endo_partitions = struct();
M_.param_names = 'h';
M_.param_names_tex = 'h';
M_.param_names_long = 'h';
M_.param_names = char(M_.param_names, 'sigma_c');
M_.param_names_tex = char(M_.param_names_tex, 'sigma\_c');
M_.param_names_long = char(M_.param_names_long, 'sigma_c');
M_.param_names = char(M_.param_names, 'beta');
M_.param_names_tex = char(M_.param_names_tex, 'beta');
M_.param_names_long = char(M_.param_names_long, 'beta');
M_.param_names = char(M_.param_names, 'adj');
M_.param_names_tex = char(M_.param_names_tex, 'adj');
M_.param_names_long = char(M_.param_names_long, 'adj');
M_.param_names = char(M_.param_names, 'tau');
M_.param_names_tex = char(M_.param_names_tex, 'tau');
M_.param_names_long = char(M_.param_names_long, 'tau');
M_.param_names = char(M_.param_names, 'gamma_p');
M_.param_names_tex = char(M_.param_names_tex, 'gamma\_p');
M_.param_names_long = char(M_.param_names_long, 'gamma_p');
M_.param_names = char(M_.param_names, 'xi_p');
M_.param_names_tex = char(M_.param_names_tex, 'xi\_p');
M_.param_names_long = char(M_.param_names_long, 'xi_p');
M_.param_names = char(M_.param_names, 'alpha');
M_.param_names_tex = char(M_.param_names_tex, 'alpha');
M_.param_names_long = char(M_.param_names_long, 'alpha');
M_.param_names = char(M_.param_names, 'gamma_w');
M_.param_names_tex = char(M_.param_names_tex, 'gamma\_w');
M_.param_names_long = char(M_.param_names_long, 'gamma_w');
M_.param_names = char(M_.param_names, 'xi_w');
M_.param_names_tex = char(M_.param_names_tex, 'xi\_w');
M_.param_names_long = char(M_.param_names_long, 'xi_w');
M_.param_names = char(M_.param_names, 'sigma_L');
M_.param_names_tex = char(M_.param_names_tex, 'sigma\_L');
M_.param_names_long = char(M_.param_names_long, 'sigma_L');
M_.param_names = char(M_.param_names, 'lambda_w');
M_.param_names_tex = char(M_.param_names_tex, 'lambda\_w');
M_.param_names_long = char(M_.param_names_long, 'lambda_w');
M_.param_names = char(M_.param_names, 'psi');
M_.param_names_tex = char(M_.param_names_tex, 'psi');
M_.param_names_long = char(M_.param_names_long, 'psi');
M_.param_names = char(M_.param_names, 'phi');
M_.param_names_tex = char(M_.param_names_tex, 'phi');
M_.param_names_long = char(M_.param_names_long, 'phi');
M_.param_names = char(M_.param_names, 'k');
M_.param_names_tex = char(M_.param_names_tex, 'k');
M_.param_names_long = char(M_.param_names_long, 'k');
M_.param_names = char(M_.param_names, 'g');
M_.param_names_tex = char(M_.param_names_tex, 'g');
M_.param_names_long = char(M_.param_names_long, 'g');
M_.param_names = char(M_.param_names, 'rho');
M_.param_names_tex = char(M_.param_names_tex, 'rho');
M_.param_names_long = char(M_.param_names_long, 'rho');
M_.param_names = char(M_.param_names, 'r_pi');
M_.param_names_tex = char(M_.param_names_tex, 'r\_pi');
M_.param_names_long = char(M_.param_names_long, 'r_pi');
M_.param_names = char(M_.param_names, 'r_Y');
M_.param_names_tex = char(M_.param_names_tex, 'r\_Y');
M_.param_names_long = char(M_.param_names_long, 'r_Y');
M_.param_names = char(M_.param_names, 'r_dpi');
M_.param_names_tex = char(M_.param_names_tex, 'r\_dpi');
M_.param_names_long = char(M_.param_names_long, 'r_dpi');
M_.param_names = char(M_.param_names, 'r_dY');
M_.param_names_tex = char(M_.param_names_tex, 'r\_dY');
M_.param_names_long = char(M_.param_names_long, 'r_dY');
M_.param_names = char(M_.param_names, 'r_k');
M_.param_names_tex = char(M_.param_names_tex, 'r\_k');
M_.param_names_long = char(M_.param_names_long, 'r_k');
M_.param_names = char(M_.param_names, 'rho_a');
M_.param_names_tex = char(M_.param_names_tex, 'rho\_a');
M_.param_names_long = char(M_.param_names_long, 'rho_a');
M_.param_names = char(M_.param_names, 'rho_I');
M_.param_names_tex = char(M_.param_names_tex, 'rho\_I');
M_.param_names_long = char(M_.param_names_long, 'rho_I');
M_.param_names = char(M_.param_names, 'rho_b');
M_.param_names_tex = char(M_.param_names_tex, 'rho\_b');
M_.param_names_long = char(M_.param_names_long, 'rho_b');
M_.param_names = char(M_.param_names, 'rho_L');
M_.param_names_tex = char(M_.param_names_tex, 'rho\_L');
M_.param_names_long = char(M_.param_names_long, 'rho_L');
M_.param_names = char(M_.param_names, 'rho_G');
M_.param_names_tex = char(M_.param_names_tex, 'rho\_G');
M_.param_names_long = char(M_.param_names_long, 'rho_G');
M_.param_names = char(M_.param_names, 'rho_pi');
M_.param_names_tex = char(M_.param_names_tex, 'rho\_pi');
M_.param_names_long = char(M_.param_names_long, 'rho_pi');
M_.param_names = char(M_.param_names, 'sd_Q');
M_.param_names_tex = char(M_.param_names_tex, 'sd\_Q');
M_.param_names_long = char(M_.param_names_long, 'sd_Q');
M_.param_names = char(M_.param_names, 'sd_p');
M_.param_names_tex = char(M_.param_names_tex, 'sd\_p');
M_.param_names_long = char(M_.param_names_long, 'sd_p');
M_.param_names = char(M_.param_names, 'sd_w');
M_.param_names_tex = char(M_.param_names_tex, 'sd\_w');
M_.param_names_long = char(M_.param_names_long, 'sd_w');
M_.param_names = char(M_.param_names, 'sd_R');
M_.param_names_tex = char(M_.param_names_tex, 'sd\_R');
M_.param_names_long = char(M_.param_names_long, 'sd_R');
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 10;
M_.endo_nbr = 27;
M_.param_nbr = 32;
M_.orig_endo_nbr = 27;
M_.aux_vars = [];
options_.varobs = cell(1);
options_.varobs(1)  = {'C'};
options_.varobs(2)  = {'I'};
options_.varobs(3)  = {'L'};
options_.varobs_id = [ 13 12 16  ];
M_.Sigma_e = zeros(10, 10);
M_.Correlation_matrix = eye(10, 10);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = 1;
M_.det_shocks = [];
options_.linear = 1;
options_.block=0;
options_.bytecode=0;
options_.use_dll=0;
M_.hessian_eq_zero = 1;
erase_compiled_function('SmetsWouters_static');
erase_compiled_function('SmetsWouters_dynamic');
M_.orig_eq_nbr = 27;
M_.eq_nbr = 27;
M_.ramsey_eq_nbr = 0;
M_.set_auxiliary_variables = exist(['./' M_.fname '_set_auxiliary_variables.m'], 'file') == 2;
M_.lead_lag_incidence = [
 1 18 0;
 2 19 45;
 3 20 46;
 4 21 0;
 5 22 0;
 6 23 0;
 0 24 0;
 7 25 47;
 8 26 48;
 9 27 0;
 0 28 49;
 10 29 50;
 11 30 51;
 12 31 0;
 0 32 52;
 0 33 0;
 13 34 0;
 0 35 53;
 0 36 0;
 14 37 0;
 0 38 54;
 15 39 55;
 16 40 56;
 0 41 0;
 0 42 57;
 0 43 0;
 17 44 0;]';
M_.nstatic = 5;
M_.nfwrd   = 5;
M_.npred   = 9;
M_.nboth   = 8;
M_.nsfwrd   = 13;
M_.nspred   = 17;
M_.ndynamic   = 22;
M_.equations_tags = {
};
M_.static_and_dynamic_models_differ = 0;
M_.exo_names_orig_ord = [1:10];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(27, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(10, 1);
M_.params = NaN(32, 1);
M_.NNZDerivatives = [125; -1; -1];
M_.params( 3 ) = 0.99;
beta = M_.params( 3 );
M_.params( 5 ) = 0.025;
tau = M_.params( 5 );
M_.params( 8 ) = 0.3;
alpha = M_.params( 8 );
M_.params( 13 ) = 5.917159763313609;
psi = M_.params( 13 );
M_.params( 6 ) = 0.469;
gamma_p = M_.params( 6 );
M_.params( 9 ) = 0.763;
gamma_w = M_.params( 9 );
M_.params( 12 ) = 0.5;
lambda_w = M_.params( 12 );
M_.params( 7 ) = 0.908;
xi_p = M_.params( 7 );
M_.params( 10 ) = 0.737;
xi_w = M_.params( 10 );
M_.params( 11 ) = 2.4;
sigma_L = M_.params( 11 );
M_.params( 2 ) = 1.353;
sigma_c = M_.params( 2 );
M_.params( 1 ) = 0.573;
h = M_.params( 1 );
M_.params( 14 ) = 1.408;
phi = M_.params( 14 );
M_.params( 4 ) = 0.1476886722788362;
adj = M_.params( 4 );
M_.params( 22 ) = 1/M_.params(3)-1+M_.params(5);
r_k = M_.params( 22 );
M_.params( 15 ) = 8.8;
k = M_.params( 15 );
M_.params( 16 ) = 0.18;
g = M_.params( 16 );
M_.params( 20 ) = 0.14;
r_dpi = M_.params( 20 );
M_.params( 19 ) = 0.099;
r_Y = M_.params( 19 );
M_.params( 21 ) = 0.159;
r_dY = M_.params( 21 );
M_.params( 17 ) = 0.961;
rho = M_.params( 17 );
M_.params( 18 ) = 1.684;
r_pi = M_.params( 18 );
M_.params( 26 ) = 0.889;
rho_L = M_.params( 26 );
M_.params( 23 ) = 0.823;
rho_a = M_.params( 23 );
M_.params( 25 ) = 0.855;
rho_b = M_.params( 25 );
M_.params( 27 ) = 0.949;
rho_G = M_.params( 27 );
M_.params( 28 ) = 0.924;
rho_pi = M_.params( 28 );
M_.params( 24 ) = 0.927;
rho_I = M_.params( 24 );
M_.params( 32 ) = 0.081;
sd_R = M_.params( 32 );
M_.params( 30 ) = 0.16;
sd_p = M_.params( 30 );
M_.params( 31 ) = 0.289;
sd_w = M_.params( 31 );
M_.params( 29 ) = 0.604;
sd_Q = M_.params( 29 );
%
% SHOCKS instructions
%
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = (1)^2;
M_.Sigma_e(2, 2) = (1)^2;
M_.Sigma_e(3, 3) = (1)^2;
M_.Sigma_e(4, 4) = (1)^2;
M_.Sigma_e(5, 5) = (1)^2;
M_.Sigma_e(6, 6) = (1)^2;
M_.Sigma_e(7, 7) = (M_.params(30))^2;
M_.Sigma_e(8, 8) = (M_.params(29))^2;
M_.Sigma_e(9, 9) = (M_.params(32))^2;
M_.Sigma_e(10, 10) = (M_.params(31))^2;
%
% INITVAL instructions
%
options_.initval_file = 0;
oo_.steady_state( 1 ) = 0;
oo_.steady_state( 2 ) = 0;
oo_.steady_state( 3 ) = 0;
oo_.steady_state( 4 ) = 0;
oo_.steady_state( 5 ) = 0;
oo_.steady_state( 6 ) = 0;
if M_.exo_nbr > 0
	oo_.exo_simul = ones(M_.maximum_lag,1)*oo_.exo_steady_state';
end
if M_.exo_det_nbr > 0
	oo_.exo_det_simul = ones(M_.maximum_lag,1)*oo_.exo_det_steady_state';
end
steady;
oo_.dr.eigval = check(M_,options_,oo_);
options_.irf = 0;
options_.order = 1;
options_.periods = 1000;
options_.simul_replic = 100;
var_list_ = char();
info = stoch_simul(var_list_);
ipred = M_.nstatic+(1:M_.nspred)';
[A,B] = kalman_transition_matrix(oo_.dr,ipred,1:M_.nspred,M_.exo_nbr);
obs_var=oo_.dr.inv_order_var(options_.varobs_id);
[C,D] = kalman_transition_matrix(oo_.dr,obs_var,1:M_.nspred,M_.exo_nbr);
temp = oo_.dr
save('SmetsWouters_results.mat', 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save('SmetsWouters_results.mat', 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save('SmetsWouters_results.mat', 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save('SmetsWouters_results.mat', 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save('SmetsWouters_results.mat', 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save('SmetsWouters_results.mat', 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save('SmetsWouters_results.mat', 'oo_recursive_', '-append');
end


disp(['Total computing time : ' dynsec2hms(toc(tic0)) ]);
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
diary off
