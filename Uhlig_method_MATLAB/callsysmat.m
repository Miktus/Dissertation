%--- Set default options and run makesysmat.m
%--- Written by Martin Schneider, Oesterreichische Nationalbank, 2007 (martin.schneider@oenb.at)
%--- Version 1.0: July 2007

if exist('SYSMAT_FIELDSIZE')~=1,
    SYSMAT_FIELDSIZE = [12 45];
end;

if exist('SYSMAT_ALIGN')~=1,
    SYSMAT_ALIGN = 'c';     % 'c': Center elements of system matrices ('l', 'r')
end;

if exist('SYSMAT_OUTPUT')~=1,
                            % 1: Writte full matrices (e.g. 'AA = [-1/rho 0 0];')
   SYSMAT_OUTPUT = 2;       % 2: Write assignment statements (e.g. 'AA(1,1) = -1/rho;')
end;

if exist('SYSMAT_DISPDIAG')~=1,
    SYSMAT_DISPDIAG = 1;    % 1: Display model diagnostics and suggest equation and variable types
end;


if exist('EQUATIONS') == 1,
    
    %--- Create system matrices
    [VARNAMES, SYSMAT_FNAM] = makesysmat(EQUATIONS, VARIABLES, SYSMAT_FIELDSIZE, SYSMAT_ALIGN, SYSMAT_OUTPUT, SYSMAT_DISPDIAG);

    %--- Refresh function and file system path cache
    rehash  
         
    %--- Evaluate system matrices (AA-NN)
    eval(SYSMAT_FNAM); 
    
end
