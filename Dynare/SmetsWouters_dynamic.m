function [residual, g1, g2, g3] = SmetsWouters_dynamic(y, x, params, steady_state, it_)
%
% Status : Computes dynamic model for Dynare
%
% Inputs :
%   y         [#dynamic variables by 1] double    vector of endogenous variables in the order stored
%                                                 in M_.lead_lag_incidence; see the Manual
%   x         [nperiods by M_.exo_nbr] double     matrix of exogenous variables (in declaration order)
%                                                 for all simulation periods
%   steady_state  [M_.endo_nbr by 1] double       vector of steady state values
%   params    [M_.param_nbr by 1] double          vector of parameter values in declaration order
%   it_       scalar double                       time period for exogenous variables for which to evaluate the model
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the dynamic model equations in order of 
%                                          declaration of the equations.
%                                          Dynare may prepend auxiliary equations, see M_.aux_vars
%   g1        [M_.endo_nbr by #dynamic variables] double    Jacobian matrix of the dynamic model equations;
%                                                           rows: equations in order of declaration
%                                                           columns: variables in order stored in M_.lead_lag_incidence followed by the ones in M_.exo_names
%   g2        [M_.endo_nbr by (#dynamic variables)^2] double   Hessian matrix of the dynamic model equations;
%                                                              rows: equations in order of declaration
%                                                              columns: variables in order stored in M_.lead_lag_incidence followed by the ones in M_.exo_names
%   g3        [M_.endo_nbr by (#dynamic variables)^3] double   Third order derivative matrix of the dynamic model equations;
%                                                              rows: equations in order of declaration
%                                                              columns: variables in order stored in M_.lead_lag_incidence followed by the ones in M_.exo_names
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

%
% Model equations
%

residual = zeros(27, 1);
T21 = (1-params(1))/((1+params(1))*params(2));
T36 = 1/(1+params(3));
T39 = params(3)/(1+params(3));
T60 = (1-params(5))/(1-params(5)+params(22));
T93 = (1-params(3)*params(7))*(1-params(7))/((1+params(3)*params(6))*params(7));
T138 = T36*(1-params(3)*params(10))*(1-params(10))/(params(10)*(1+1/params(12)*(1+params(12))*params(11)));
T142 = params(2)/(1-params(1));
lhs =y(30);
rhs =params(1)/(1+params(1))*y(11)+1/(1+params(1))*y(51)-T21*(y(31)-y(47))+T21*(y(20)-y(46));
residual(1)= lhs-rhs;
lhs =y(29);
rhs =T36*y(10)+T39*y(50)+params(4)/(1+params(3))*y(28)+(params(3)*y(45)-y(19))/(1+params(3));
residual(2)= lhs-rhs;
lhs =y(28);
rhs =(-(y(31)-y(47)))+T60*y(49)+params(22)/(1-params(5)+params(22))*y(52)+x(it_, 8);
residual(3)= lhs-rhs;
lhs =y(27);
rhs =(1-params(5))*y(9)+y(10)*params(5);
residual(4)= lhs-rhs;
lhs =y(25);
rhs =y(47)*params(3)/(1+params(3)*params(6))+params(6)/(1+params(3)*params(6))*y(7)+T93*(params(8)*y(32)+(1-params(8))*y(26)-y(18)+x(it_, 7));
residual(5)= lhs-rhs;
lhs =y(26);
rhs =T39*y(48)+T36*y(8)+y(47)*T39-y(25)*(1+params(3)*params(9))/(1+params(3))+y(7)*params(9)/(1+params(3))-T138*(y(26)-params(11)*y(33)-T142*(y(30)-params(1)*y(11))-y(21)-x(it_, 10));
residual(6)= lhs-rhs;
lhs =y(33);
rhs =y(9)+(-y(26))+y(32)*(1+params(13));
residual(7)= lhs-rhs;
lhs =y(34);
rhs =y(18)*params(14)+y(9)*params(8)*params(14)+y(32)*params(13)*params(8)*params(14)+y(33)*(1-params(8))*params(14);
residual(8)= lhs-rhs;
lhs =y(34);
rhs =y(30)*(1-params(5)*params(15)-params(16))+y(29)*params(5)*params(15)+y(22);
residual(9)= lhs-rhs;
lhs =y(31);
rhs =params(17)*y(12)+(1-params(17))*(y(23)+params(18)*(y(7)-y(23))+params(19)*(y(34)-y(44)))+params(20)*(y(25)-y(7))+params(21)*(y(34)-y(44)-(y(13)-y(17)))+x(it_, 9);
residual(10)= lhs-rhs;
lhs =y(40);
rhs =T21*(y(20)-y(46))+params(1)/(1+params(1))*y(16)+1/(1+params(1))*y(56)-T21*(y(41)-y(53));
residual(11)= lhs-rhs;
lhs =y(39);
rhs =(params(3)*y(45)-y(19))/(1+params(3))+T36*y(15)+T39*y(55)+params(4)/(1+params(3))*y(38);
residual(12)= lhs-rhs;
lhs =y(38);
rhs =(-(y(41)-y(53)))+T60*y(54)+params(22)/(1-params(5)+params(22))*y(57);
residual(13)= lhs-rhs;
lhs =y(37);
rhs =(1-params(5))*y(14)+params(5)*y(15);
residual(14)= lhs-rhs;
residual(15) = y(35);
lhs =0;
rhs =params(8)*y(42)+(1-params(8))*y(36)-y(18);
residual(16)= lhs-rhs;
lhs =y(36);
rhs =params(11)*y(43)+T142*(y(40)-params(1)*y(16))-y(21);
residual(17)= lhs-rhs;
lhs =y(43);
rhs =y(14)+(-y(36))+(1+params(13))*y(42);
residual(18)= lhs-rhs;
lhs =y(44);
rhs =y(18)*params(14)+params(8)*params(14)*y(14)+params(13)*params(8)*params(14)*y(42)+(1-params(8))*params(14)*y(43);
residual(19)= lhs-rhs;
lhs =y(44);
rhs =y(22)+(1-params(5)*params(15)-params(16))*y(40)+params(5)*params(15)*y(39);
residual(20)= lhs-rhs;
lhs =y(19);
rhs =params(24)*y(2)+x(it_, 2);
residual(21)= lhs-rhs;
lhs =y(20);
rhs =params(25)*y(3)+x(it_, 3);
residual(22)= lhs-rhs;
lhs =y(21);
rhs =params(26)*y(4)+x(it_, 4);
residual(23)= lhs-rhs;
lhs =y(22);
rhs =params(27)*y(5)+x(it_, 5);
residual(24)= lhs-rhs;
lhs =y(18);
rhs =params(23)*y(1)+x(it_, 1);
residual(25)= lhs-rhs;
lhs =y(23);
rhs =params(28)*y(6)+x(it_, 6);
residual(26)= lhs-rhs;
residual(27) = y(24);
if nargout >= 2,
  g1 = zeros(27, 67);

  %
  % Jacobian matrix
  %

  g1(1,20)=(-T21);
  g1(1,46)=T21;
  g1(1,47)=(-T21);
  g1(1,11)=(-(params(1)/(1+params(1))));
  g1(1,30)=1;
  g1(1,51)=(-(1/(1+params(1))));
  g1(1,31)=T21;
  g1(2,19)=(-((-1)/(1+params(3))));
  g1(2,45)=(-T39);
  g1(2,28)=(-(params(4)/(1+params(3))));
  g1(2,10)=(-T36);
  g1(2,29)=1;
  g1(2,50)=(-T39);
  g1(3,47)=(-1);
  g1(3,28)=1;
  g1(3,49)=(-T60);
  g1(3,31)=1;
  g1(3,52)=(-(params(22)/(1-params(5)+params(22))));
  g1(3,65)=(-1);
  g1(4,9)=(-(1-params(5)));
  g1(4,27)=1;
  g1(4,10)=(-params(5));
  g1(5,18)=T93;
  g1(5,7)=(-(params(6)/(1+params(3)*params(6))));
  g1(5,25)=1;
  g1(5,47)=(-(params(3)/(1+params(3)*params(6))));
  g1(5,26)=(-(T93*(1-params(8))));
  g1(5,32)=(-(T93*params(8)));
  g1(5,64)=(-T93);
  g1(6,21)=(-T138);
  g1(6,7)=(-(params(9)/(1+params(3))));
  g1(6,25)=(1+params(3)*params(9))/(1+params(3));
  g1(6,47)=(-T39);
  g1(6,8)=(-T36);
  g1(6,26)=1-(-T138);
  g1(6,48)=(-T39);
  g1(6,11)=T138*(-(T142*(-params(1))));
  g1(6,30)=T138*(-T142);
  g1(6,33)=T138*(-params(11));
  g1(6,67)=(-T138);
  g1(7,26)=1;
  g1(7,9)=(-1);
  g1(7,32)=(-(1+params(13)));
  g1(7,33)=1;
  g1(8,18)=(-params(14));
  g1(8,9)=(-(params(8)*params(14)));
  g1(8,32)=(-(params(13)*params(8)*params(14)));
  g1(8,33)=(-((1-params(8))*params(14)));
  g1(8,34)=1;
  g1(9,22)=(-1);
  g1(9,29)=(-(params(5)*params(15)));
  g1(9,30)=(-(1-params(5)*params(15)-params(16)));
  g1(9,34)=1;
  g1(10,23)=(-((1-params(17))*(1-params(18))));
  g1(10,7)=(-((1-params(17))*params(18)-params(20)));
  g1(10,25)=(-params(20));
  g1(10,12)=(-params(17));
  g1(10,31)=1;
  g1(10,13)=params(21);
  g1(10,34)=(-(params(21)+(1-params(17))*params(19)));
  g1(10,17)=(-params(21));
  g1(10,44)=(-((1-params(17))*(-params(19))-params(21)));
  g1(10,66)=(-1);
  g1(11,20)=(-T21);
  g1(11,46)=T21;
  g1(11,53)=(-T21);
  g1(11,16)=(-(params(1)/(1+params(1))));
  g1(11,40)=1;
  g1(11,56)=(-(1/(1+params(1))));
  g1(11,41)=T21;
  g1(12,19)=(-((-1)/(1+params(3))));
  g1(12,45)=(-T39);
  g1(12,38)=(-(params(4)/(1+params(3))));
  g1(12,15)=(-T36);
  g1(12,39)=1;
  g1(12,55)=(-T39);
  g1(13,53)=(-1);
  g1(13,38)=1;
  g1(13,54)=(-T60);
  g1(13,41)=1;
  g1(13,57)=(-(params(22)/(1-params(5)+params(22))));
  g1(14,14)=(-(1-params(5)));
  g1(14,37)=1;
  g1(14,15)=(-params(5));
  g1(15,35)=1;
  g1(16,18)=1;
  g1(16,36)=(-(1-params(8)));
  g1(16,42)=(-params(8));
  g1(17,21)=1;
  g1(17,36)=1;
  g1(17,16)=(-(T142*(-params(1))));
  g1(17,40)=(-T142);
  g1(17,43)=(-params(11));
  g1(18,36)=1;
  g1(18,14)=(-1);
  g1(18,42)=(-(1+params(13)));
  g1(18,43)=1;
  g1(19,18)=(-params(14));
  g1(19,14)=(-(params(8)*params(14)));
  g1(19,42)=(-(params(13)*params(8)*params(14)));
  g1(19,43)=(-((1-params(8))*params(14)));
  g1(19,44)=1;
  g1(20,22)=(-1);
  g1(20,39)=(-(params(5)*params(15)));
  g1(20,40)=(-(1-params(5)*params(15)-params(16)));
  g1(20,44)=1;
  g1(21,2)=(-params(24));
  g1(21,19)=1;
  g1(21,59)=(-1);
  g1(22,3)=(-params(25));
  g1(22,20)=1;
  g1(22,60)=(-1);
  g1(23,4)=(-params(26));
  g1(23,21)=1;
  g1(23,61)=(-1);
  g1(24,5)=(-params(27));
  g1(24,22)=1;
  g1(24,62)=(-1);
  g1(25,1)=(-params(23));
  g1(25,18)=1;
  g1(25,58)=(-1);
  g1(26,6)=(-params(28));
  g1(26,23)=1;
  g1(26,63)=(-1);
  g1(27,24)=1;

if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],27,4489);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],27,300763);
end
end
end
end
