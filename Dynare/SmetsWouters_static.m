function [residual, g1, g2, g3] = SmetsWouters_static(y, x, params)
%
% Status : Computes static model for Dynare
%
% Inputs : 
%   y         [M_.endo_nbr by 1] double    vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1] double     vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1] double   vector of parameter values in declaration order
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the static model equations 
%                                          in order of declaration of the equations.
%                                          Dynare may prepend or append auxiliary equations, see M_.aux_vars
%   g1        [M_.endo_nbr by M_.endo_nbr] double    Jacobian matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%   g2        [M_.endo_nbr by (M_.endo_nbr)^2] double   Hessian matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%   g3        [M_.endo_nbr by (M_.endo_nbr)^3] double   Third derivatives matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

residual = zeros( 27, 1);

%
% Model equations
%

T19 = (1-params(1))/((1+params(1))*params(2));
T30 = 1/(1+params(3));
T32 = params(3)/(1+params(3));
T51 = (1-params(5))/(1-params(5)+params(22));
T80 = (1-params(3)*params(7))*(1-params(7))/((1+params(3)*params(6))*params(7));
T122 = T30*(1-params(3)*params(10))*(1-params(10))/(params(10)*(1+1/params(12)*(1+params(12))*params(11)));
lhs =y(13);
rhs =y(13)*params(1)/(1+params(1))+y(13)*1/(1+params(1))-T19*(y(14)-y(8));
residual(1)= lhs-rhs;
lhs =y(12);
rhs =y(12)*T30+y(12)*T32+params(4)/(1+params(3))*y(11)+(params(3)*y(2)-y(2))/(1+params(3));
residual(2)= lhs-rhs;
lhs =y(11);
rhs =(-(y(14)-y(8)))+y(11)*T51+params(22)/(1-params(5)+params(22))*y(15)+x(8);
residual(3)= lhs-rhs;
lhs =y(10);
rhs =(1-params(5))*y(10)+y(12)*params(5);
residual(4)= lhs-rhs;
lhs =y(8);
rhs =y(8)*params(3)/(1+params(3)*params(6))+y(8)*params(6)/(1+params(3)*params(6))+T80*(y(15)*params(8)+(1-params(8))*y(9)-y(1)+x(7));
residual(5)= lhs-rhs;
lhs =y(9);
rhs =T32*y(9)+T30*y(9)+y(8)*T32-y(8)*(1+params(3)*params(9))/(1+params(3))+y(8)*params(9)/(1+params(3))-T122*(y(9)-params(11)*y(16)-params(2)/(1-params(1))*(y(13)-y(13)*params(1))-y(4)-x(10));
residual(6)= lhs-rhs;
lhs =y(16);
rhs =y(10)+(-y(9))+y(15)*(1+params(13));
residual(7)= lhs-rhs;
lhs =y(17);
rhs =y(1)*params(14)+y(10)*params(8)*params(14)+y(15)*params(13)*params(8)*params(14)+y(16)*(1-params(8))*params(14);
residual(8)= lhs-rhs;
lhs =y(17);
rhs =y(13)*(1-params(5)*params(15)-params(16))+y(12)*params(5)*params(15)+y(5);
residual(9)= lhs-rhs;
lhs =y(14);
rhs =y(14)*params(17)+(1-params(17))*(y(6)+params(18)*(y(8)-y(6))+params(19)*(y(17)-y(27)))+x(9);
residual(10)= lhs-rhs;
lhs =y(23);
rhs =params(1)/(1+params(1))*y(23)+1/(1+params(1))*y(23)-T19*(y(24)-y(18));
residual(11)= lhs-rhs;
lhs =y(22);
rhs =(params(3)*y(2)-y(2))/(1+params(3))+T30*y(22)+T32*y(22)+params(4)/(1+params(3))*y(21);
residual(12)= lhs-rhs;
lhs =y(21);
rhs =(-(y(24)-y(18)))+T51*y(21)+params(22)/(1-params(5)+params(22))*y(25);
residual(13)= lhs-rhs;
lhs =y(20);
rhs =(1-params(5))*y(20)+params(5)*y(22);
residual(14)= lhs-rhs;
residual(15) = y(18);
lhs =0;
rhs =params(8)*y(25)+(1-params(8))*y(19)-y(1);
residual(16)= lhs-rhs;
lhs =y(19);
rhs =params(11)*y(26)+params(2)/(1-params(1))*(y(23)-params(1)*y(23))-y(4);
residual(17)= lhs-rhs;
lhs =y(26);
rhs =y(20)+(-y(19))+(1+params(13))*y(25);
residual(18)= lhs-rhs;
lhs =y(27);
rhs =y(1)*params(14)+params(8)*params(14)*y(20)+params(13)*params(8)*params(14)*y(25)+(1-params(8))*params(14)*y(26);
residual(19)= lhs-rhs;
lhs =y(27);
rhs =y(5)+(1-params(5)*params(15)-params(16))*y(23)+params(5)*params(15)*y(22);
residual(20)= lhs-rhs;
lhs =y(2);
rhs =y(2)*params(24)+x(2);
residual(21)= lhs-rhs;
lhs =y(3);
rhs =y(3)*params(25)+x(3);
residual(22)= lhs-rhs;
lhs =y(4);
rhs =y(4)*params(26)+x(4);
residual(23)= lhs-rhs;
lhs =y(5);
rhs =y(5)*params(27)+x(5);
residual(24)= lhs-rhs;
lhs =y(1);
rhs =y(1)*params(23)+x(1);
residual(25)= lhs-rhs;
lhs =y(6);
rhs =y(6)*params(28)+x(6);
residual(26)= lhs-rhs;
residual(27) = y(7);
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(27, 27);

  %
  % Jacobian matrix
  %

T321 = 1-(params(1)/(1+params(1))+1/(1+params(1)));
  g1(1,8)=(-T19);
  g1(1,13)=T321;
  g1(1,14)=T19;
  g1(2,2)=(-((params(3)-1)/(1+params(3))));
  g1(2,11)=(-(params(4)/(1+params(3))));
  g1(2,12)=1-(T30+T32);
  g1(3,8)=(-1);
  g1(3,11)=1-T51;
  g1(3,14)=1;
  g1(3,15)=(-(params(22)/(1-params(5)+params(22))));
  g1(4,10)=1-(1-params(5));
  g1(4,12)=(-params(5));
  g1(5,1)=T80;
  g1(5,8)=1-(params(3)/(1+params(3)*params(6))+params(6)/(1+params(3)*params(6)));
  g1(5,9)=(-(T80*(1-params(8))));
  g1(5,15)=(-(T80*params(8)));
  g1(6,4)=(-T122);
  g1(6,8)=(-(params(9)/(1+params(3))+T32-(1+params(3)*params(9))/(1+params(3))));
  g1(6,9)=1-(T30+T32-T122);
  g1(6,13)=T122*(-((1-params(1))*params(2)/(1-params(1))));
  g1(6,16)=T122*(-params(11));
  g1(7,9)=1;
  g1(7,10)=(-1);
  g1(7,15)=(-(1+params(13)));
  g1(7,16)=1;
  g1(8,1)=(-params(14));
  g1(8,10)=(-(params(8)*params(14)));
  g1(8,15)=(-(params(13)*params(8)*params(14)));
  g1(8,16)=(-((1-params(8))*params(14)));
  g1(8,17)=1;
  g1(9,5)=(-1);
  g1(9,12)=(-(params(5)*params(15)));
  g1(9,13)=(-(1-params(5)*params(15)-params(16)));
  g1(9,17)=1;
  g1(10,6)=(-((1-params(17))*(1-params(18))));
  g1(10,8)=(-((1-params(17))*params(18)));
  g1(10,14)=1-params(17);
  g1(10,17)=(-((1-params(17))*params(19)));
  g1(10,27)=(-((1-params(17))*(-params(19))));
  g1(11,18)=(-T19);
  g1(11,23)=T321;
  g1(11,24)=T19;
  g1(12,2)=(-((params(3)-1)/(1+params(3))));
  g1(12,21)=(-(params(4)/(1+params(3))));
  g1(12,22)=1-(T30+T32);
  g1(13,18)=(-1);
  g1(13,21)=1-T51;
  g1(13,24)=1;
  g1(13,25)=(-(params(22)/(1-params(5)+params(22))));
  g1(14,20)=1-(1-params(5));
  g1(14,22)=(-params(5));
  g1(15,18)=1;
  g1(16,1)=1;
  g1(16,19)=(-(1-params(8)));
  g1(16,25)=(-params(8));
  g1(17,4)=1;
  g1(17,19)=1;
  g1(17,23)=(-((1-params(1))*params(2)/(1-params(1))));
  g1(17,26)=(-params(11));
  g1(18,19)=1;
  g1(18,20)=(-1);
  g1(18,25)=(-(1+params(13)));
  g1(18,26)=1;
  g1(19,1)=(-params(14));
  g1(19,20)=(-(params(8)*params(14)));
  g1(19,25)=(-(params(13)*params(8)*params(14)));
  g1(19,26)=(-((1-params(8))*params(14)));
  g1(19,27)=1;
  g1(20,5)=(-1);
  g1(20,22)=(-(params(5)*params(15)));
  g1(20,23)=(-(1-params(5)*params(15)-params(16)));
  g1(20,27)=1;
  g1(21,2)=1-params(24);
  g1(22,3)=1-params(25);
  g1(23,4)=1-params(26);
  g1(24,5)=1-params(27);
  g1(25,1)=1-params(23);
  g1(26,6)=1-params(28);
  g1(27,7)=1;
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],27,729);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],27,19683);
end
end
end
end
