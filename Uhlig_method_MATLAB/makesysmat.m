function [VARNAMES, SYSMAT_FNAM] = makesysmat(EQUATIONS, VARIABLES, SYSMAT_FIELDSIZE, SYSMAT_ALIGN, SYSMAT_OUTPUT, SYSMAT_DISPDIAG)
% makesysmat: Generates coefficient matrices AA-NN for Uhlig's toolkit from model equations
%
% This function creates MATLAB code which defines the system matrices AA - NN 
% from the model equations and writes it into file SYSMAT_FNAM
%--------------------------------------------------------------------------
% % Written by: Martin Schneider, Oesterreichische Nationalbank, 2007
%               martin.schneider@oenb.at
%--------------------------------------------------------------------------
% This code can be used, copied and modified by everybody at his/her own
% risk, but not sold or otherwise used for commercial purposes.
%--------------------------------------------------------------------------
% Version 1.0: July 2007
% This version was written with MATLAB 7.0.4 and tested with MATLAB versions 7.0.1 and 6.5
% There is no warranty that the programme is compatible with other MATLAB versions
%--------------------------------------------------------------------------
%  Usage: Copy the files callsysmat.m and makesysmat.m in the toolkit folder
%         and add 'callsysmat' as the first command in the file 'do_it.m'
%--------------------------------------------------------------------------     
%  Inputs:
%  EQUATIONS         [N,3]: Equation type ~ Equation name ~ equation (cell array)
%                           Equation type: 'D' ... Deterministic equation
%                                          'E' ... Expectational equation 
%                                          'S' ... Shock equation
%                                          ''  ... Skip equation
%  VARIABLES         [M,3]: Variable name (case sensitive) ~ Variable type ~ Description (cell array)
%                           Variable type: 'X' ... Endogenous state variable 
%                                          'Y' ... Endogenous other (jump) variable  
%                                          'Z' ... Exogenous state variable 
%                                          'U' ... Innovation to exogenous state variable
%                                          ''  ... Skip variable
%--------------------------------------------------------------------------     
%  Options:
%  SYSMAT_FIELDSIZE  [1,1]: Size of matrix elements 
%                           (can also be 2x1 containing SYSMAT_FIELDSIZE for 
%                            SYSMAT_OUTPUT = 1 and SYSMAT_OUTPUT = 2)
%  SYSMAT_ALIGN      [Str]: Alignment of matrix elements  ('l', 'r', 'c')
%  SYSMAT_OUTPUT     [1x1]: 1 ... Write full matrices
%                           2 ... Write assignment statements (e.g. 'AA(1,1) = -1/rho;')     
%  SYSMAT_DISPDIAG   [1x1]: 1 ... Display model diagnostics and suggest equation
%                                 and variable types
%--------------------------------------------------------------------------     
%  Output:   
%  VARNAMES       : Variable names as character array 
%  SYSMAT_FNAM    : Name of file in which system matrices are written
%--------------------------------------------------------------------------     
% Structure of the model in Uhlig's notation:
%   0      = AA x(t) + BB x(t-1) + CC y(t) + DD z(t)
%   0      = E_t [ FF x(t+1) + GG x(t) + HH x(t-1) + JJ y(t+1) + KK y(t) + LL z(t+1) + MM z(t)]
%   z(t+1) = NN z(t) + epsilon(t+1),    with E_t [ epsilon(t+1) ] = 0
%
%   x(t)      : Endogenous state variables 
%   y(t)      : Endogenous other variables
%   z(t)      : Exogenous state variables
%   epsilon(t): Innovation to exogenous state variable
%
% Syntax of the equations:
%   - Examples for valid time subscripts are x(t-1), x(-1), x(0), x(t),
%     x(t+1), x(1), x(+1)
%   - Omit the expectation operator , i.e. 
%         0 = E_t [-eta*c(t+1)+ r(t+1)+eta*c(t)] ==>
%         0 = -eta*c(t+1)+ r(t+1)+eta*c(t)
%   - You can omit the '=' sign. Then all variables are placed on the RHS of the equation
%   - Multiple occurences of variables per equation are not allowed
%   - You may include valid functions in your expression (MATLAB functions
%     or functions included in the search path).
%   - Shock equations must have time subscript t+1 on the LHS and t on the
%     RHS for the exogenous variables
%
%--------------------------------------------------------------------------
% Example:
%
% 
%     VARIABLES = {...
%                 'k'         'X'  'Capital'  
%                 'c'         'Y'  'Consumption'
%                 'r'         'Y'  'Interest rate'
%                 'y'         'Y'  'Output'
%                 'z'         'Z'  'Technology shock'
%                 'epsilon'   'U'  'Innovation to shock'};
%                
% 
%     EQUATIONS  = {...
%     'D' '1: Capital costs' ...
%             '0 = - (1 - betta *(1-delta))*(1-rho)* k(t-1) + (1 - betta*(1-delta))*z(t) - r(t)'
%     
%     'D' '2: Evolution of capital stock' ...
%             '0 = - K_bar/C_bar*k(t) + K_bar/(betta*C_bar)*k(t-1) + (1 + delta*K_bar/C_bar)*z(t) - c(t)'
% 
%     'D' '3: Production function' ...
%             '0 = - y(t) + z(t) + rho * k(t-1)'
%     
%     'E' '4: Consumption euler equation' ...
%             '0 = - eta*c(t+1) + r(t+1) + eta*c(t)'
%     
%     'S' '5: Evolution of technology shock' ...
%             'z(t+1) =   psi*z(t) + epsilon(t+1)'};
%
%
%--------------------------------------------------------------------------------------------------------
% ==> Writes the following code in file SYSMAT_MODELNAME.m (if SYSMAT_OUTPUT == 1)
% 
% %--- DETERMINISTIC EQUATIONS 
% 
% %---            k(t)    
% AA = [           0         % 1: Capital costs
%             -K_bar/C_bar    % 2: Evolution of capital stock
%                  0      ]; % 3: Production function
% 
% %---           k(t-1)   
% BB = [      -(1-betta*(1-delta))*(1-rho)    % 1: Capital costs
%             K_bar/(betta*C_bar)    % 2: Evolution of capital stock
%                 rho     ]; % 3: Production function
% 
% %---            c(t)        r(t)        y(t)    
% CC = [           0          -1           0         % 1: Capital costs
%                 -1           0           0         % 2: Evolution of capital stock
%                  0           0          -1      ]; % 3: Production function
% 
% %---            z(t)    
% DD = [      (1-betta*(1-delta))    % 1: Capital costs
%             (1+delta*K_bar/C_bar)    % 2: Evolution of capital stock
%                  1      ]; % 3: Production function
% 
% 
% %--- EXPECTATIONAL EQUATIONS 
% 
% %---           k(t+1)   
% FF = [           0      ]; % 4: Consumption euler equation
% 
% %---            k(t)    
% GG = [           0      ]; % 4: Consumption euler equation
% 
% %---           k(t-1)   
% HH = [           0      ]; % 4: Consumption euler equation
% 
% %---           c(t+1)      r(t+1)      y(t+1)   
% JJ = [         -eta          1           0      ]; % 4: Consumption euler equation
% 
% %---            c(t)        r(t)        y(t)    
% KK = [          eta          0           0      ]; % 4: Consumption euler equation
% 
% %---           z(t+1)   
% LL = [           0      ]; % 4: Consumption euler equation
% 
% %---            z(t)    
% MM = [           0      ]; % 4: Consumption euler equation
% 
% 
% %--- SHOCK PROCESSES 
% 
% %---        epsilon(t+1)
% NN = [          psi     ]; % 5: Evolution of technology shock
% 
%
%--------------------------------------------------------------------------


%--- Extract variable names for figures (pad with blanks)
% VARNAMES = VARIABLES(:,3); 
VARNAMES = [];
maxlen = 0;
for i=1:size(VARIABLES,1)  % Determine longest variable name
    if length(VARIABLES{i,3}) > maxlen
        maxlen = length(VARIABLES{i,3});
    end
end

for i=1:size(VARIABLES,1)
    if strcmp(VARIABLES{i,2},'') == 0
        name = VARIABLES{i,3};
        diff = maxlen - length(name);
        if diff > 0
            name = [name char(zeros(1,diff)+32)];
        end
        VARNAMES = [VARNAMES; name];
    end
end

%--- Extract information from input matrices

    EqDetIdx = [];
    EqExpIdx = [];
    EqShkIdx = [];
    for r = 1:size(EQUATIONS,1)
        eqtyp = upper(EQUATIONS{r,1});
        if strcmp(    eqtyp,'D') == 1
            EqDetIdx = [EqDetIdx; r];
        elseif strcmp(eqtyp,'E') == 1
            EqExpIdx = [EqExpIdx; r];
        elseif strcmp(eqtyp,'S') == 1
            EqShkIdx = [EqShkIdx; r];
        elseif strcmp(eqtyp,'') == 1
        else
            error(['Error in equation matrix for equation ''' upper(EQUATIONS{r,2}) ''' : Invalid equation type ''' eqtyp ''' (use ''D(eterministic)'', ''E(xpectational)'', ''S(hock)'') or '''' (to omit equation)'])
        end        
    end


    xNames   = {};
    yNames   = {};
    zNames   = {};
    shkNames = {};

    for r = 1:size(VARIABLES,1)
        varnam = VARIABLES{r,1};
        vartyp = upper(VARIABLES{r,2});
        if strcmp(vartyp,'X') == 1
            xNames   = [xNames; {varnam}];
        elseif strcmp(vartyp,'Y') == 1
            yNames   = [yNames; {varnam}];
        elseif strcmp(vartyp,'Z') == 1
            zNames   = [zNames; {varnam}];
        elseif strcmp(vartyp,'U') == 1
            shkNames = [shkNames; {varnam}];
        elseif strcmp(vartyp,'') == 1
        else
            error(['Error in variable matrix for variable ''' varnam ''': Invalid variable type ''' vartyp ''' (use ''X(Endog. state variable)'', ''Y(jump variable)'', ''Z(exog. state variable)'') or ''U'' (shock innovation)'])
        end        
    end


    EQTYPE  = EQUATIONS(:,1);
    EqNames = EQUATIONS(:,2);
    EQUATIONS   = EQUATIONS(:,3);


    if length(SYSMAT_FIELDSIZE) == 2
        SYSMAT_FIELDSIZE = SYSMAT_FIELDSIZE(SYSMAT_OUTPUT);
    end


    VarNames = [xNames; yNames; zNames];


    EqDet      = EQUATIONS(EqDetIdx);
    EqExp      = EQUATIONS(EqExpIdx);
    EqShk      = EQUATIONS(EqShkIdx);
    
    EqDetNames = EqNames(EqDetIdx);    
    EqExpNames = EqNames(EqExpIdx);
    EqShkNames = EqNames(EqShkIdx);
    
    NrEqDet   = length(EqDetIdx);
    NrEqExp   = length(EqExpIdx);
    NrEqShk   = length(EqShkIdx);
    
    NrVarState     = size(xNames,1);
    NrVarJump = size(yNames,1);
    NrVarExo  = size(zNames,1);

    crlf = char([13 10]);

    %--- Initialize matrices
    if SYSMAT_OUTPUT == 1
        
        SysMat.AA = num2cell(zeros(NrEqDet, NrVarState));
        SysMat.BB = num2cell(zeros(NrEqDet, NrVarState));
        SysMat.CC = num2cell(zeros(NrEqDet, NrVarJump));
        SysMat.DD = num2cell(zeros(NrEqDet, NrVarExo));

        SysMat.FF = num2cell(zeros(NrEqExp, NrVarState));
        SysMat.GG = num2cell(zeros(NrEqExp, NrVarState));
        SysMat.HH = num2cell(zeros(NrEqExp, NrVarState));

        SysMat.JJ = num2cell(zeros(NrEqExp, NrVarJump));
        SysMat.KK = num2cell(zeros(NrEqExp, NrVarJump));
        SysMat.LL = num2cell(zeros(NrEqExp, NrVarExo));
        SysMat.MM = num2cell(zeros(NrEqExp, NrVarExo));
        SysMat.NN = num2cell(zeros(NrEqShk, NrEqShk));
        
    else
        
        str = '';
        str = [str 'AA = zeros(' int2str(NrEqDet) ', ' int2str(NrVarState) ');' crlf];
        str = [str 'BB = zeros(' int2str(NrEqDet) ', ' int2str(NrVarState) ');' crlf];
        str = [str 'CC = zeros(' int2str(NrEqDet) ', ' int2str(NrVarJump) ');'  crlf];
        str = [str 'DD = zeros(' int2str(NrEqDet) ', ' int2str(NrVarExo) ');' crlf];
        
        str = [str 'FF = zeros(' int2str(NrEqExp) ', ' int2str(NrVarState) ');' crlf];
        str = [str 'GG = zeros(' int2str(NrEqExp) ', ' int2str(NrVarState) ');' crlf];
        str = [str 'HH = zeros(' int2str(NrEqExp) ', ' int2str(NrVarState) ');' crlf];

        str = [str 'JJ = zeros(' int2str(NrEqExp) ', ' int2str(NrVarJump) ');' crlf];
        str = [str 'KK = zeros(' int2str(NrEqExp) ', ' int2str(NrVarJump) ');' crlf];
        str = [str 'LL = zeros(' int2str(NrEqExp) ', ' int2str(NrVarExo) ');' crlf];
        str = [str 'MM = zeros(' int2str(NrEqExp) ', ' int2str(NrVarExo) ');' crlf];
        str = [str 'NN = zeros(' int2str(NrEqShk) ', ' int2str(NrEqShk) ');' crlf];
    end

    eqwritten = zeros(length(EQUATIONS),1);
   
    %----------------------------------------------------------------------
    %--- Big loop over all equations
    %----------------------------------------------------------------------
    okall        =  1;
    ErrorTextAll = '';
    diagmat      = {};
    for EqIdx = 1: size(EQUATIONS, 1)    
    
        %--- Convert equation from string to cell array
        EqCell = {};
        rem = strrep(EQUATIONS{EqIdx},' ', '');
        
        %--- Index of equation within deterministic / expectational / shock equations
        EqType = lower(EQTYPE{EqIdx,1});
        if strcmp(EqType, 'd') == 1
            EqTypeLong = 'deterministic';
            EqNr       = find(EqDetIdx==EqIdx);
        elseif strcmp(EqType, 'e') == 1
            EqTypeLong = 'expectational';
            EqNr       = find(EqExpIdx==EqIdx);
        elseif strcmp(EqType, 's') == 1
            EqTypeLong = 'shock';
            EqNr       = find(EqShkIdx==EqIdx);
        else
            %--- Skip equation
            EqNr = 0;
        end
        EqNr = EqNr(1);
        
        
        if strcmp(EqType,'') == 0
            
            eqwritten(EqIdx) = 1;            
            EqTextError = ['Error in ' EqTypeLong ' equation ' EqNames{EqIdx} ':' char(13) '"' strrep(EQUATIONS{EqIdx},' ','') '"'];
            
            %--- Parse equations for proper syntax
            [ok, errstr, diagmat] = parser(rem, EqIdx, EqType, xNames, yNames, zNames, shkNames, diagmat);

            if not(ok)
                %--- Concatenate error message
                ErrorTextAll = [ErrorTextAll EqTextError char(13) errstr  char(13) char(13)];
                okall = 0;
            end

            if ok
                %--- Split equations into variable names and coefficients
                [coeffs, varnames] = spliteq(rem, EqType, xNames, yNames, zNames, shkNames, EqTextError);

                if SYSMAT_OUTPUT == 1
                    %--- Write coefficients into system matrices
                    for i = 1:length(varnames)
                        SysMat = writeinallmat(varnames{i}, coeffs{i}, EqType, EqNr, xNames, yNames, zNames, shkNames, SysMat);
                    end
                else
                    %--- Create string with assignment statements
                    str = writeassignment(varnames, coeffs, EqType, EqNr, EqNames{EqIdx}, xNames, yNames, zNames, SYSMAT_FIELDSIZE, str);

                end
            end
        end % if  strcmp(EqType,'') ...
            
    end % loop over equations

    if SYSMAT_DISPDIAG
        %-- Display diagnostics           
        showdiagnostics(diagmat, VARIABLES, EQUATIONS, EqNames, EQTYPE, NrEqDet, NrEqExp, NrEqShk, NrVarState, NrVarJump, NrVarExo);
    else
        disp(' ')
        disp('You can display model diagnostics by setting SYSMAT_DISPDIAG=1')
        disp(' ')
    end
    
    if okall == 0;
        %--- Show errors and diagnostic matrices

        
        %--- Show errors and diagnostic matrices
        error(ErrorTextAll)
        return
    end
    
    %----------------------------------------------------------------------
    %--- Generate string
    %----------------------------------------------------------------------

    if SYSMAT_OUTPUT == 1
        %--- Write out whole matrices
        str = [    crlf '%--- DETERMINISTIC EQUATIONS ', crlf crlf];

        str = [str writemat('AA', SysMat.AA, xNames , EqDetNames,  0, SYSMAT_FIELDSIZE, SYSMAT_ALIGN)];       
        str = [str writemat('BB', SysMat.BB, xNames , EqDetNames, -1, SYSMAT_FIELDSIZE, SYSMAT_ALIGN)];
        str = [str writemat('CC', SysMat.CC, yNames , EqDetNames,  0, SYSMAT_FIELDSIZE, SYSMAT_ALIGN)];
        str = [str writemat('DD', SysMat.DD, zNames , EqDetNames,  0, SYSMAT_FIELDSIZE, SYSMAT_ALIGN)];

        str = [str crlf '%--- EXPECTATIONAL EQUATIONS ', crlf crlf];
        str = [str writemat('FF', SysMat.FF, xNames , EqExpNames,  1, SYSMAT_FIELDSIZE, SYSMAT_ALIGN)];
        str = [str writemat('GG', SysMat.GG, xNames , EqExpNames,  0, SYSMAT_FIELDSIZE, SYSMAT_ALIGN)];
        str = [str writemat('HH', SysMat.HH, xNames , EqExpNames, -1, SYSMAT_FIELDSIZE, SYSMAT_ALIGN)];
        str = [str writemat('JJ', SysMat.JJ, yNames , EqExpNames,  1, SYSMAT_FIELDSIZE, SYSMAT_ALIGN)];
        str = [str writemat('KK', SysMat.KK, yNames , EqExpNames,  0, SYSMAT_FIELDSIZE, SYSMAT_ALIGN)];
        str = [str writemat('LL', SysMat.LL, zNames , EqExpNames,  1, SYSMAT_FIELDSIZE, SYSMAT_ALIGN)];
        str = [str writemat('MM', SysMat.MM, zNames , EqExpNames,  0, SYSMAT_FIELDSIZE, SYSMAT_ALIGN)];

        str = [str crlf '%--- SHOCK PROCESSES ', crlf crlf];
        str = [str writemat('NN', SysMat.NN, shkNames , EqShkNames,  1, SYSMAT_FIELDSIZE, SYSMAT_ALIGN)];
    end

    %----------------------------------------------------------------------
    %--- Write into file
    %----------------------------------------------------------------------
    
    %--- Get name of model file
    [ST,I]      = dbstack;
    SYSMAT_FNAM = ST(end).name;
    
    %--- delete '.m' from file name (for versions prior to 7.0) 
    idx = strfind2(lower(SYSMAT_FNAM),'.m');
    if isempty(idx) == 0
        SYSMAT_FNAM = SYSMAT_FNAM(1:idx-1);
    end
        
    %--- Separate file name from path (for versions prior to 7.0)
    idx = max(find(SYSMAT_FNAM =='\'));
    if isempty(idx) == 1
        %--- Look for '/'
        idx = max(find(SYSMAT_FNAM =='/'));
    end

    if isempty(idx) == 0
        %--- SYSMAT_FNAM contains full path --> separate
        SYSMAT_FNAM = SYSMAT_FNAM(idx+1:end);
    end

    
    SYSMAT_FNAM = [SYSMAT_FNAM  '_SYSMAT'];

    fid = fopen([SYSMAT_FNAM '.m'], 'w');
    fprintf(fid, '%s', str);     
    fclose(fid);

    %----------------------------------------------------------------------  
    %--- Display message
    %----------------------------------------------------------------------  
    
    disp(['The following equations were written into coefficient matrices in file "' SYSMAT_FNAM '"'])
    for i = 1:length(eqwritten)            
        if eqwritten(i)
            disp(EqNames{i})
        end
    end
    
    disp(char(13))
    
    if sum(1-eqwritten) > 0
        disp('Warning: you have omitted the following equation(s):')
        for i = 1:length(eqwritten)            
            if not(eqwritten(i))
                disp(EqNames{i})
            end
        end
    else
        disp('All equations written!')
    end
    
    

function [ok, errstr, diagmat] = parser(eq, EqIdx, EqType, xNames, yNames, zNames, shkNames, diagmat)
%**************************************************************************
% Parse formula and check for syntax errors
%**************************************************************************

    vorg        = '(';
    ok         = 1;
    errstr     = '';
    BrackLevel = 0;
    InFunction = 0;
    BrackLevelFunction = 0;
    oplist1    = '+-*/^=';
    oplist2    = '+-*/^,=';
    oplist     = oplist1;
    varlist    = {};
    eqlst      = {};
    eq         = strrep(eq, ' ', '');
    
    %--- Check for multiple '=' signs
    idx = strfind2(eq ,'=');
    if length(idx) > 1
        ok = 0;
        errstr = 'multiple "=" signs'
        return
    end
    
    
    %-- Generate list of all valid variable names for that equation
    allnames = validnames(EqType, xNames, yNames, zNames, shkNames);
    
    i = 1;

    while i <= length(eq)
        
        x = eq(i);
        
        if strcmp(vorg,'#function#')
            %------------------------------------------------------
            %--- After MATLAB function: operand or '(' expected ............ to be done
            %------------------------------------------------------
           
            if strcmp('(',x)
                BrackLevel = BrackLevel + 1;
                vorg = '(';
            else
                ok = 0;
                errstr =  ['"(" after MATLAB function "' op '" expected'];
                return                     
            end
            
        elseif not(isempty(strfind2('()', vorg))) | not(isempty(strfind2([oplist '()'], vorg)))            
            %------------------------------------------------------
            %--- After operator or '(': operand or '(' expected
            %------------------------------------------------------
            
            if strcmp('(',vorg) & strcmp('-',x)
                
                vorg = '-';
                
            elseif strcmp('(',x)
                
                vorg = '(';
                BrackLevel = BrackLevel + 1;
                if InFunction
                    BrackLevelFunction = BrackLevelFunction + 1;
                end
                
            elseif strcmp(')',x)
                
                vorg = ')';
                BrackLevel = BrackLevel - 1;

                if BrackLevel < 0;
                    ok = 0;
                    errstr = '"(" missing ';
                    return
                end
                
                if InFunction
                    BrackLevelFunction = BrackLevelFunction - 1;
                    if BrackLevel == 0
                        InFunction = 0;
                        oplist = oplist1;
                    end
                end
 
                
            elseif isvarname(x) |  not(isempty(str2num(x)))
                
            %------------------------------------------------------
            %--- Check operand
            %------------------------------------------------------

                op = x;
                if isvarname(op)
                    
                    %-- Check if expression is valid name for MATLAB or model variable
                    while isvarname(op) & (i < length(eq))
                        i  = i + 1;
                        x  = eq(i);
                        op = [op x];
                    end
                    
                    if not(isvarname(op))
                        i = i - 1;
                        op = op(1:end-1);
                    end

                    if isfunction(op) 
                        %--- Function
                        vorg = '#function#';
                        InFunction = 1;
                        oplist = oplist2;
                        
                    else
                        %--- Check for valid variable
                        i2 = i + 1;
                        if i2 <= length(eq)

                            if strcmp(eq(i2),'(')

                                BrackLevel = BrackLevel + 1;

                                %--- Eventually variable: Look for time subscript
                                while not(strcmp(eq(i2),')')) & i2 < length(eq)
                                   op = [op eq(i2)];
                                   i2 = i2 + 1;
                                end

                                if strcmp(eq(i2),')')
                                    BrackLevel = BrackLevel - 1;
                                    if InFunction
                                        BrackLevelFunction = BrackLevelFunction - 1;
                                    end
                                end
                                op = [op eq(i2)];

                                %--- Append to diagnostic matrix
                                % = Variable ~ lead/lag ~ Equation nr.
                                %--- Separate varnam from lag
                                ixx    = strfind2(op, '(');
                                varnam = op(1:ixx-1);
                                lag    = str2num(strrep(op(ixx+1:end-1),'t',''));
                                if isempty(lag)
                                    lag = 0; 
                                end
                                diagmat = [diagmat; {varnam lag EqIdx op}];

                                %--- Check for valid lag of variable
                                if isvariable(op, EqType, xNames, yNames, zNames, shkNames)
                                    %--- check for multiple occurences of variables
                                    

                                    if indcv(op, varlist)
                                        ok = 0; 
                                        errstr = ['multiple occurence of variable ""' op '"'];
                                        return 
                                    else
                                        varlist = [varlist; {op}];
                                        
                                    end
                                    i = i2-0; 
                                else
                                    
                                    ok = 0; 
                                    errstr = ['Invalid variable (or invalid lead/lag for this equation type) "' op '"'];
                                    return 
                                end % if isvariable(op ...                                                                      

                            end % if strcmp(eq(i2),'(')
                        end
                        x    = op;
                        vorg = '#operand#';   

                    end % if isfunction(SYSMAT_FNAM)

                                            
                elseif not(isempty(str2num(op))) 
                    
                    %-- Number
                    while not(isempty(str2num(op))) & (i < length(eq))
                        i  = i + 1;
                        x  = eq(i);
                        op = [op x];
                    end
                    
                    if isempty(str2num(op))
                        i = i - 1;
                        op = op(1:end-1);
                    end
                    
                    x    = op;
                    vorg = '#operand#';   
                    
                end % if isvarname(op)

            elseif not(strfind2('+-',x))
                
                ok = 0; 
                errstr = ['Operator or ")" expected after "' eq(1:i-1) '"'];
                return                
                
            end
            
        elseif strcmp(vorg, ')') | strcmp(vorg, '#operand#')
            
            %--- After operand or '(': operator or ')' or time subscript e.g '(t+1)'expected 
            
            if not(isempty(strfind2(oplist, x)))
                %--- Operator
                vorg = x;
            elseif strcmp(x,')')
                vorg = ')';
                BrackLevel = BrackLevel - 1;

                if BrackLevel < 0;
                    ok = 0;
                    errstr = '"(" missing ';
                    return
                end
                
                if InFunction
                    BrackLevelFunction = BrackLevelFunction - 1;
                    if BrackLevel == 0
                        InFunction = 0;
                        oplist = oplist1;
                    end
                end
                 
            else
                if strcmp(x,'(')
                    %--- Check, if expression is a variable name (e.g 'x(t)')
                    i2 = strfind2(eq(i+1:end),')');
                    if isempty(i2)
                        ok = 0;
                        errstr = ['Operator or ")" expected after "' eq(1:i-1) '"'];
                        return                    
                    end
                    if not(isvariable(op, EqType, xNames, yNames, zNames, shkNames))
                        ok = 0;
                        errstr = ['Unknown variable or invalid lag "' op '"'];
                        return                    
                    end
                    i = i + i2(1);
                else
                    ok = 0;
                    errstr = ['Operator or ")" expected after "' eq(1:i-1) '"'];
                    return                    
                end
                
            end
        else
            ok = 0;
            errstr = ['Invalid character "' x '"'];
            return                                
        end % if 
        
        eqlst = [eqlst; {x}];    % append to list 
        i = i + 1;
    end % while

    if BrackLevel > 0
        ok = 0;
        errstr = '")" missing';
        return  
    elseif not(isempty(strfind2(oplist, vorg)))
        ok = 0;
        errstr = 'Incomplete expression';
        return  
    else
    end
    

function [coeffs, varnames, errtxt] = spliteq(equation, EqType, xNames, yNames, zNames, shkNames, EqTextError)
%**************************************************************************
% Split up equation in coefficients and varnames
% e.g. '0 = x(t) + (1-a)*(A*y(t)-z(t)' -->
%       coefficients   varnames
%            1           x(t)
%          (1-a)*A       y(t)
%         -(1-a)         r(t)
%**************************************************************************

    coeffs    = {};
    varnames  = {};
    errtxtall = '';
    
    equation = strrep(equation, ' ', '');

    
    [addterms, errtxt] = splitaddterm(equation);
    
    if not(isempty(errtxt))        
        errtxt = [EqTextError char(13) errtxt];
        return
    end

    
    for i = 1:length(addterms)
        [coeffs, varnames, vartermret, errtxt] = multiplyterm(addterms{i}, coeffs, varnames, EqType, xNames, yNames, zNames, shkNames);
        if not(isempty(errtxt))
        errtxt = [EqTextError char(13) errtxt];
        end
    end
    

function [coeffs, varnames, vartermret, errtxt] = multiplyterm(term, coeffs, varnames, EqType, xNames, yNames, zNames, shkNames)
%**************************************************************************
% Multiply out terms to get terms that are linear in the variables
%
% e.g. (1+A)*(y(t)+B*z(t))*(1-delta)
%            coeff            varnam
% -->   (1+A)*(1-delta)        y(t)
%       (1+a)*B*(1-delta)      z(t)  
%
%**************************************************************************

    errtxt     = '';
    vartermret = '';
    varterm    = '';

    %-- Generate list of all valid variable names for that equation
    allnames = validnames(EqType, xNames, yNames, zNames, shkNames);

    %--- Split multiplicative terms    
    [multterms, errtxt] = splitmultterm(term, allnames);


    if isempty(multterms)
         if not(strcmp(term,'0') | strcmp(term,'-0'))
            errtxt = ['term "' term '" contains no variable or variable has invalid lag'];
            return
        end
    end
    if not(isempty(errtxt))
        return
    end
 
    isvar = zeros(size(multterms,1),1);


    for i = 1:size(multterms,1)

        VNamAll  = 0;
        
        if not(containsvar(multterms{i}, allnames))
            %--- No variable contained: term includes only coefficients'
            isvar(i) = 0;
        else
            %--- Variable contained
            isvar(i) = 1;
            
            %--- Split each multiplicative term in additive terms
            [addterms, errtxt] = splitaddterm(multterms{i});
            if not(isempty(errtxt))
                return
            end
        
                        
            %--- Recursive call: For each additive term multiply out
            
            varterm = '';

            for j = 1:size(addterms,1)
                
                if not(containsvar(addterms{j}, allnames))
                    errtxt = ['Additive term "' addterms{j} '" does not contain a variable and would be lost.' ]
                    return
                end
                
                if not(isatom(addterms{j}, allnames))
                
                    [coeffs, varnames, vartermret, errtxt] = multiplyterm(addterms{j}, ...
                     coeffs, varnames, EqType, xNames, yNames, zNames, shkNames);

                    if not(isempty(errtxt))
                        return
                    end

                    if not(isempty(vartermret))
                        if strcmp(vartermret(1),'-')
                            varterm = [varterm vartermret];
                        else
                            varterm = [varterm '+' vartermret];
                        end
                    end
                else
                    %--- addterm{j} is variable --> append
                    adt     = addterms{j};
                    if strcmp(adt(1),'-')
                        varterm = [varterm addterms{j}];
                    else
                        varterm = [varterm '+' addterms{j}];
                    end
                end % if not isatom
            end % for j = 1:
            
            if strcmp(varterm(1),'+')
                varterm = varterm(2:end);
            end
            multterms{i} = varterm;

        end % if not containsvar(...                   
    end % for i = 1:size(multterms,1)

    
    if sum(isvar) == 0       
         vartermret = '';
         return
    elseif sum(isvar) > 1
        errtxt = ['term "' term '" contains a non-linear combination of variables'];
        return
    end
        

    %--- Multiply out coefficients
    c = '';
    for i = 1:size(multterms,1) 
        if not(isvar(i))
            trm = multterms{i};
            if strcmp(trm(1),'/')
                c = [c trm];
            else
                c = [c '*' trm];
            end
        end
    end
    c = c(2:end);
    if isempty(c)
        c = '1';
    end

    %--- Generate term that contains only variable names of additive terms !!!
    idx     = find(isvar);
    idx     = idx(1);
    varterm = multterms{idx};
    
    %--- Write mutiplied coefficients for all variables in 'coeff'
    
    %-- Split up additive term with variables
    [splittedvarterms, errtxt] = splitaddterm(varterm);
    if not(isempty(errtxt))
        return
    end
    
     %--- Multiply (all) variable(s) in (additive) term with coefficients
    vartermret = '';
    for i = 1:size(splittedvarterms,1)
        
        a = splittedvarterms{i};
       
        if length(a) > 1 & strcmp(a(1),'-')
            %--- variable with '-': Flig sign of coeff
            vartermret = [vartermret a];
            
            a = a(2:end);
            if strcmp(c(1),'-')
                c2 = c(2:end);
            else
                c2 = ['-' c];
            end
        else
            c2 = c;
            vartermret = [vartermret '+' a];
        end


        %--- Search for variable in varnames
        idx = indcv(a,varnames);
        if idx == 0
            %--- Not found --> append row
            varnames = [varnames; {a} ];
            coeffs   = [coeffs  ; {c2}];
        else
            %--- Already existent --> Multiply

            cold = coeffs{idx};
            if strcmp(cold,'1')
                coeffs{idx} = [c2];
            elseif strcmp(cold,'-1')
                coeffs{idx} = ['-' c2];
                
            elseif strcmp(cold(1),'-')                
                %--- Flip sign
                if strcmp(c2(1),'-')
                    coeffs{idx} = [c2(2:end) '*' cold(2:end)];
                else
                    coeffs{idx} = ['-' c2 '*' cold(2:end)];
                end
            else
                coeffs{idx} = [c2 '*' cold];                
            end
        end
    end
       
    if strcmp(vartermret(1),'+')
        vartermret = vartermret(2:end);
    end
    
    
function ok = isatom(term, varnames)
%**************************************************************************
% Checks, whether multiplicative term is atomar
%**************************************************************************
 
    ok = 1;
    
    if strcmp(term(1),'-')
        term = term(2:end);
    end

    if indcv(term, varnames)
        %--- valid variable name
        return
    end
    
    for i = 1:length(term)
        if strfind2('+-*/^',term(i))
            ok = 0;
            return
        end
    end


function [multterms, errtxt] = splitmultterm(term, varlst)
%**************************************************************************
% Split multiplicative term into components
%**************************************************************************


    multterms = {};
    errtxt    = '';
    
    if size(term) == 0
        return
    end
    
    if not(containsvar(term,varlst))
        %--- Return, if no variable is included
        return
    end

    token     = '';
    brackets  = 0;
    term      = strrep(term,' ', '');
    denum     = 0; % denumerator (after /)
    
    i = 1;
    for i = 1:length(term)

        x = term(i);
        token = [token x];

        if strfind2('({[',x) > 0            
            brackets = brackets + 1;
        elseif strfind2(')}]',x) > 0            
            brackets = brackets - 1;
            if brackets == -1
                errtxt = ['"(" missing in term "' term '"'];
                return
            end
        elseif strfind2('*',x) > 0
            if brackets == 0
                if length(token) > 1
                    if denum
                        if containsvar(token(1:end-1),varlst)
                            errtxt = ['variable in denominator in expression "' term '"'];
                            return
                        else
                            if size(multterms,1) == 0
                                multterms = [multterms ; {['1/' token(1:end-1)]}];
                            else
                                multterms = [multterms ; {['/'  token(1:end-1)]}];
                            end
                        end
                    else
                        if not(strcmp(token(1:end-1),'1'))
                            multterms = [multterms ; {token(1:end-1)}];
                        end
                    end
                    token     = '';
                else
                   errtxt = ['error in term "' term '"'];
                   return
                end
                denum = 0;
            end
            
        elseif strfind2('/',x) > 0
            if brackets == 0
                if length(token) > 1
                    if denum
                        if containsvar(token(1:end-1),varlst)
                            errtxt = ['variable in denominator in expression "' term '"'];
                            return
                        else
                            if size(multterms,1) == 0
                                multterms = [multterms ; {['1/' token(1:end-1)]}];
                            else
                                multterms = [multterms ; {['/' token(1:end-1)]}];
                            end
                        end
                    else
                        if not(strcmp(token(1:end-1),'1'))
                            multterms = [multterms ; {token(1:end-1)}];
                        end
                    end
                    token     = '';
                else                   
                   errtxt = ['error in term "' term '"'];
                   return
                end
                denum = 1;
            end
        end
        i = i + 1;
    end
    
    if length(token) > 0
        if denum
            if containsvar(token(1:end),varlst)
                errtxt = ['variable in denominator in expression "' term '"'];
                return
            else
                if size(multterms,1) == 0
                    multterms = [multterms ; {['1/' token(1:end)]}];
                else
                    multterms = [multterms ; {['/' token(1:end)]}];
                end
            end
        else
            if not(strcmp(token(1:end),'1'))
                multterms = [multterms ; {token(1:end)}];
            end
        end
    end
    

    if brackets > 0
        errtxt =  ['")" missing in term "' term '"'];
        return
    end

    %--- Remove leading and trailing brackets
    if size(multterms,1) == 1  &  strcmp(term(1),'(') & strcmp(term(end),')')
        term = term(2:end-1);
        if checkbrackets(term)
            [multterms, errtxt] = splitmultterm(term, varlst);
        end
    end

    
function [addterms, errtxt] = splitaddterm(term)
%**************************************************************************
% Split additive term into components
%**************************************************************************


    addterms = {};
    errtxt   = '';
    
    if size(term) == 0
        return
    end

    idx = strfind2(term ,'=');
    if not(isempty(idx))
        LHS = 1;
    else
        LHS = 0;
    end
    
    token     = '';
    brackets  = 0;
    term      = strrep(term,' ', '');
     
    i = 1;
    for i = 1:length(term)

        x = term(i);
 
        if strfind2('({[',x) > 0
            
            token = [token x];
            brackets = brackets + 1;
            
        elseif strfind2(')}]',x) > 0
            
            token = [token x];
            brackets = brackets - 1;
            if brackets == -1
                addterms = {};
                errtxt =  ['"(" missing in term "' term '"'];
                return
            end
            
        elseif strfind2('+-=',x) > 0 | (i == length(term))
            
             if (brackets == 0) & not(strcmp(token,''))
                %--- Term finished: Append to list

                if LHS
                    %--- Flip sign of current term
                    if length(token) > 1 
                        if strcmp(token(1),'-')
                            token = token(2:end);
                        elseif strcmp(token(1),'+')
                            token = ['-' token(2:end)];
                        else
                            token = ['-' token];
                        end
                    else
                        token = ['-' token];
                    end
                else
                    %--- RHS: Trim leading '+'
                    if length(token) > 1 & strcmp(token(1),'+')
                        token = token(2:end);
                    end
                end
 
                %--- Append token to addterms
                if not(strcmp(token,'0'))
                    addterms = [addterms ; {token}];
                end
                
                %--- Initialize next term
                if strcmp(x,'=')
                    LHS = 0;
                    token = '';
                end
                
                %--- Sign for next term
                if strcmp(x,'-')
                    token  = '-';
                else
                    token = '';
                end
                
             else
                 %--- continue current term
                 token = [token x];
             end %  if (brackets == ...
             
        else
            %--- continue current term
            token = [token x];
        end        
        i = i + 1;
    end
                   
    if i == length(term) + 1
        %--- Append last token                
        if not(strcmp(token,'0'))
            addterms = [addterms ; {token}];
        end
    end

    if brackets > 0
        errtxt = (['")" expected in term "' term '"']);
        return
    end

    %--- Remove leading and trailing brackets
    if size(addterms,1) == 1  &  strcmp(term(1),'(') & strcmp(term(end),')')
        term = term(2:end-1);
        %--- Check for brackets
        if checkbrackets(term)
            [addterms, errtxt] = splitaddterm(term);
        end
    end

    
function ok = containsvar(expr, allnames)
%**************************************************************************
% Checks, whether string or cell array (vector) contains variable
%**************************************************************************

   ok = 0;

   if isstr(expr)
        %--- String
        for j = 1:length(allnames)
            if strfind2(expr, allnames{j})
                ok = 1;
                return
            end
        end
    else
        %--- Character array
        for i = 1:size(expr,1)
            for j = 1:length(allnames)
                if strfind2(expr{i}, allnames{j})
                    ok = 1;
                    return
                end
            end        
        end
    end


function ok = checkbrackets(term)
%**************************************************************************
% Check expression for brackets
%**************************************************************************
    
    ok = 1;
    lvl = 0;
    for i = 1:length(term)
        if strcmp(term(i),'(')
            lvl = lvl + 1;
        elseif strcmp(term(i),')')
            if lvl == 0 
                ok = 0;
                break;
            else
                lvl = lvl - 1;
            end
        end
    end
    if lvl > 0
        ok = 0
    end
    

function ok = isfunction(name)
%**************************************************************************
% Determine if name is a MATLAB function (built-in or user-defined)
%**************************************************************************

    if indcv(name,{'i'; 'pi' ; 'psi'}) > 0
        ok = 0;
        return
    end

    if isvarname(name)

        eval(['fh=functions(@' name ');']);

        if isempty(fh.file)
            ok = 0;
        else
            ok = 1;
        end
    else
        ok = 0;
    end
                             
function [str] = writemat(name,  mat, varnames, eqtext, lag, SYSMAT_FIELDSIZE, type)
%**************************************************************************
% Write matrix including header and equation names to string
%**************************************************************************
    
    SYSMAT_FIELDSIZE2 = SYSMAT_FIELDSIZE;

    crlf = char([13 10]);
    
    if lag > 0 
        lagtxt = ['(t+' num2str(lag) ')'];
    elseif lag < 0 
        lagtxt = ['(t' num2str(lag) ')'];
    else
        lagtxt = '(t)';
    end

    
    %--- Write variable names    
    str = [fillblanks('%---', SYSMAT_FIELDSIZE2, 'l')];
    if length(varnames) == 1
        varnames = repmat(varnames, size(mat,2),1);
    end
    for i = 1:size(varnames,1)
        str = [str fillblanks([varnames{i} lagtxt], SYSMAT_FIELDSIZE, type)] ;
    end
    str = [str crlf];


    %--- Write matrix elements
    if size(mat,1) == 0
        %--- Empty matrix
        str = [str fillblanks([name ' = [];'],SYSMAT_FIELDSIZE, 'l')];
    else
       %--- Non-empty matrix
        for r = 1:size(mat,1)
            %--- Begin of line
            if r == 1
                linestr = [fillblanks([name ' = [ '], SYSMAT_FIELDSIZE2, 'l')];
            else
                linestr = [fillblanks('', SYSMAT_FIELDSIZE2-1, type) ' '] ;
            end
            %--- Write matrix elements
            for c = 1:size(mat,2)
                x = mat{r,c};
                if x == 0
                    x = '0';
                end
                rem = mod(length(linestr),SYSMAT_FIELDSIZE);
                linestr = [linestr fillblanks(x, SYSMAT_FIELDSIZE-1-rem, type) ' '] ;            
            end
            %-- End of line
            if r == size(mat,1)
                linestr = [linestr ']; % ' eqtext{r} crlf];
            else
                linestr = [linestr '   % ' eqtext{r} crlf] ;            
            end
            str = [str linestr];

        end
    end

    str =[str crlf];


function str = fillblanks(str, SYSMAT_FIELDSIZE, type)
%**************************************************************************
% Fill up string with blanks
%**************************************************************************

    l = length(str);
    if l >= SYSMAT_FIELDSIZE
        return
    end
    
    if strcmp(type,'l')
        str = [str char(zeros(1,SYSMAT_FIELDSIZE-l)+32)];
    elseif strcmp(type,'r')
        str = [char(zeros(1,SYSMAT_FIELDSIZE-l)+32) str];
    elseif strcmp(type,'c')
        l = SYSMAT_FIELDSIZE-length(str);
        if mod(l,2) == 0
            l1 = l/2;
            l2 = l/2;
        else
            l1 = floor(l/2);
            l2 = ceil(l/2);
        end
        
        str = [char(zeros(1,l1)+32) str char(zeros(1,l2)+32)];

    end
    
    
function [ridx, cidx] = findvarnam(varnam, vnamlst, lag)
%**************************************************************************
% Search varnam with given lag in vnamlst
%**************************************************************************
    
    for ridx = 1:size(vnamlst,1)
        for cidx = 1:size(vnamlst,2)
            if lag == 0
                x1 =   vnamlst{ridx,cidx};
                x2 =  [vnamlst{ridx,cidx} '(t)'];
                x3 =  [vnamlst{ridx,cidx} '(0)'];
            elseif lag > 0
                x1 =  [vnamlst{ridx,cidx} '(+' num2str(lag) ')'];
                x2 =  [vnamlst{ridx,cidx} '(t+' num2str(lag) ')'];
                x3 =  [vnamlst{ridx,cidx} '(' num2str(lag) ')'];
            else
                x1 =  [vnamlst{ridx,cidx} '(t' num2str(lag) ')'];
                x2 =  [vnamlst{ridx,cidx} '(' num2str(lag) ')'];
                x3 = x1;
            end
            if strcmp(varnam,x1) == 1
                return
            elseif strcmp(varnam,x2) == 1
                return
            elseif strcmp(varnam,x3) == 1
                return
            end
        end
    end
    
    ridx = 0;
    cidx = 0;

    
function  validnames = validnames(EqType, xNames, yNames, zNames, shkNames)
%**************************************************************************
% make list of all valid variables
%**************************************************************************

    validnames = {};
    vnamlst    = [xNames; yNames; zNames; shkNames];
    typ        = [repmat(cellstr('x'),size(xNames,1),1)  ; ...
                  repmat(cellstr('y'),size(yNames,1),1)  ; ...
                  repmat(cellstr('z'),size(zNames,1),1)  ; ...
                  repmat(cellstr('s'),size(shkNames,1),1); ...
                 ];
    
    for idx = 1:size(vnamlst,1)
        typ{idx};
        if strcmp(EqType,'')
            validnames =  [validnames; ...
                              {[vnamlst{idx} '(t-1)']};
                              {[vnamlst{idx} '(-1)' ]}; 
                              {[vnamlst{idx} '(t)'  ]};
                              {[vnamlst{idx} '(0)'  ]};
                              {[vnamlst{idx} '(t+1)']};
                              {[vnamlst{idx} '(1)'  ]};
                              {[vnamlst{idx} '(+1)' ]}
                          ];
                      
        elseif strcmp(EqType,'d')
            
            %--- Deterministic equations
            if strcmp(typ{idx},'x')
                validnames =  [validnames; ...
                                  {[vnamlst{idx} '(t-1)']};
                                  {[vnamlst{idx} '(-1)' ]};
                                  {[vnamlst{idx} '(t)'  ]};
                                  {[vnamlst{idx} '(0)'  ]}
                              ];
            elseif strcmp(typ{idx},'y')
                validnames =  [validnames; ...
                                  {[vnamlst{idx} '(t)']};
                                  {[vnamlst{idx} '(0)']}
                              ];
            elseif strcmp(typ{idx},'z')
                validnames =  [validnames; ...
                                  {[vnamlst{idx} '(t)']};
                                  {[vnamlst{idx} '(0)']}
                              ];
            end
            
        elseif strcmp(EqType,'e')
            
            %--- Expectational equations
            
            if strcmp(typ{idx},'x')
                validnames =  [validnames; ...
                                  {[vnamlst{idx} '(t-1)']};
                                  {[vnamlst{idx} '(-1)' ]};
                                  {[vnamlst{idx} '(t)'  ]};
                                  {[vnamlst{idx} '(0)'  ]};
                                  {[vnamlst{idx} '(t+1)']};
                                  {[vnamlst{idx} '(+1)' ]};
                                  {[vnamlst{idx} '(1)'  ]}
                              ];
            elseif strcmp(typ{idx},'y')
                validnames =  [validnames; ...
                                  {[vnamlst{idx} '(t)'  ]};
                                  {[vnamlst{idx} '(0)'  ]};
                                  {[vnamlst{idx} '(t+1)']};
                                  {[vnamlst{idx} '(+1)' ]};
                                  {[vnamlst{idx} '(1)'  ]}
                              ];
            elseif strcmp(typ{idx},'z')
                validnames =  [validnames; ...
                                  {[vnamlst{idx} '(t)'  ]};
                                  {[vnamlst{idx} '(0)'  ]};
                                  {[vnamlst{idx} '(t+1)']};
                                  {[vnamlst{idx} '(+1)' ]};
                                  {[vnamlst{idx} '(1)'  ]}
                              ];
            end
            
        elseif strcmp(EqType,'s')
        
            %--- Shock equations
            if strcmp(typ{idx},'z')
                validnames =  [validnames; ...
                                  {[vnamlst{idx} '(t)'  ]};
                                  {[vnamlst{idx} '(0)'  ]};
                                  {[vnamlst{idx} '(t+1)']};
                                  {[vnamlst{idx} '(+1)' ]};
                                  {[vnamlst{idx} '(1)'  ]};
                              ];
            elseif strcmp(typ{idx},'s')
                validnames =  [validnames; ...
                                  {[vnamlst{idx} '(t+1)']};
                                  {[vnamlst{idx} '(+1)' ]};
                                  {[vnamlst{idx} '(1)'  ]}
                              ];
            end
            
        end
        
    end
    

function ok = isvariable(VarNam, EqType, xNames, yNames, zNames, shkNames)
%**************************************************************************
% Determine, if variable name is valid in the sense of Uhlig
%**************************************************************************

   
    ok = 1;

    if strcmp(EqType, 'd')
        %--- Deterministic

        %--- AA: xNames{t}
        idx = findvarnam(VarNam, xNames, 0);
        if idx > 0
            return 
        end

        %--- BB: xNames{t-1}
        idx = findvarnam(VarNam, xNames, -1);
        if idx > 0
            return
        end

        %--- CC: yNames{t}
        idx = findvarnam(VarNam, yNames, 0);
        if idx > 0
            return
        end

        %--- DD: zNames{t}
        idx = findvarnam(VarNam, zNames, 0);
        if idx > 0
            return
        end

    elseif strcmp(EqType, 'e')
        
        %--- Expectational
        
        %--- FF: xNames{t+1}
        idx = findvarnam(VarNam, xNames, 1);
        if idx > 0
            return
        end

        %--- GG: xNames{t}
        idx = findvarnam(VarNam, xNames, 0);
        if idx > 0
            return
        end

        %--- HH: xNames{t-1}
        idx = findvarnam(VarNam, xNames, -1);
        if idx > 0
            return
        end

        %--- JJ: yNames{t+1}
        idx = findvarnam(VarNam, yNames, 1);
        if idx > 0
            return
        end

        %--- KK: yNames{t}
        idx = findvarnam(VarNam, yNames, 0);
        if idx > 0
            return
        end

        %--- LL: zNames{t+1}
        idx = findvarnam(VarNam, zNames, 1);
        if idx > 0
            return
        end

        %--- MM: yNames{t}
        idx = findvarnam(VarNam, zNames, 0);
        if idx > 0
            return
        end
        
    else
        
        %--- LHS of shock processes: zNames{t+1}
        idx = findvarnam(VarNam, zNames, 1);
        if idx > 0
            return
        end

        %--- NN: zNames{t}
        idx = findvarnam(VarNam, zNames, 0);
        if idx > 0
            return
        end
        %--- Innovation in shock processes: shkNames{t+1}
        idx = findvarnam(VarNam, shkNames, 1);
        if idx > 0
            return
        end
        
    end

    ok = 0;

    
function [mat, ok] = writeintomat(mat, varnam, coeff, vnamlst, eqnum, lag)
%**************************************************************************
% Write coefficient in matrix element identified by varnam and eqnum
%**************************************************************************

    ok = 0;
    idx = findvarnam(varnam, vnamlst, lag);
    
    if idx > 0
        mat{eqnum, idx} = coeff;
        ok = 1;
    end

    
function [SysMat] = ... 
          writeinallmat(varnam, coeff, EqType, EqNr, xNames, yNames, zNames, shkNames, SysMat)
%**************************************************************************
% Write coefficient in all matrices
%**************************************************************************

    written = 0;

    if strcmp(EqType, 'd');
    %--- Deterministic equations

        %--- AA: xNames{t}
        [SysMat.AA, ok] = writeintomat(SysMat.AA, varnam, coeff, xNames, EqNr, 0);
        written = written + ok;

        %--- BB: xNames{t-1}
        [SysMat.BB, ok] = writeintomat(SysMat.BB, varnam, coeff, xNames, EqNr, -1);
        written = written + ok;

        %--- CC: yNames{t}
        [SysMat.CC, ok] = writeintomat(SysMat.CC, varnam, coeff, yNames, EqNr, 0);
        written = written + ok;

        %--- DD: zNames{t}
        [SysMat.DD, ok] = writeintomat(SysMat.DD, varnam, coeff, zNames, EqNr, 0);
        written = written + ok;

    elseif strcmp(EqType, 'e');
        %--- Expectational equations

        %--- FF: xNames{t+1}
        [SysMat.FF, ok] = writeintomat(SysMat.FF, varnam, coeff, xNames, EqNr, 1);
        written = written + ok;

        %--- GG: xNames{t}
        [SysMat.GG, ok] = writeintomat(SysMat.GG, varnam, coeff, xNames, EqNr, 0);
        written = written + ok;

        %--- HH: xNames{t-1}
        [SysMat.HH, ok] = writeintomat(SysMat.HH, varnam, coeff, xNames, EqNr, -1);
        written = written + ok;

        %--- JJ: yNames{t+1}
        [SysMat.JJ, ok] = writeintomat(SysMat.JJ, varnam, coeff, yNames, EqNr, 1);
        written = written + ok;

        %--- KK: yNames{t}
        [SysMat.KK, ok] = writeintomat(SysMat.KK, varnam, coeff, yNames, EqNr, 0);
        written = written + ok;

        %--- LL: zNames{t+1}
        [SysMat.LL, ok] = writeintomat(SysMat.LL, varnam, coeff, zNames, EqNr, 1);
        written = written + ok;

        %--- MM: zNames{t}
        [SysMat.MM, ok] = writeintomat(SysMat.MM, varnam, coeff, zNames, EqNr, 0);
        written = written + ok;

    else
        %--- NN: Autoregressive matrix for z{t}
        [SysMat.NN, ok] = writeintomat(SysMat.NN, varnam, coeff, zNames, EqNr, 0);
        
        %--- Check zNames(t+1) and shkNames(t+1) 
        if ok
            written = ok;
        else            
            idx = findvarnam(varnam, zNames, 1) + findvarnam(varnam, shkNames, 1);
            if idx == 0
                written = 0;
            else
                written = 1;
            end
        end

    end % if EqType ...
    
    if not(written)
        error(['Error during writing into system matrices: Invalid variable name "' varnam '"'])
    end

    
function str = writeassignment(varnames, coeffs, EqType, EqNr, EqNam, xNames, yNames, zNames, SYSMAT_FIELDSIZE, str);
%**************************************************************************
% Write assignment statements
%**************************************************************************

    crlf = char([13]);
    
    str = [str crlf crlf '%--- ' EqNam crlf crlf];
    
    
    for i = 1:length(varnames)
        
        matnam = '';
        
        %--- Split name and lag of variable
        vn   = varnames{i};
        idx1 = strfind2(vn, '(');
        idx2 = strfind2(vn, ')');
        lag  = vn(idx1+1:idx2-1);
        lag  = strrep(lag,'t', '');
        if isempty(lag)
            lag = 0;
        else
            lag = str2num(lag);
        end        
        vn  = vn(1:idx1-1);

        %--- Get name of matrix and index of variable
        if strcmp(EqType, 'd');         %--- Deterministic equation
            
            typ   = 'deterministic';
            
            %--- xNames
            idx = indcv(vn, xNames);
            if idx > 0
                VarNr = idx;
                if lag == 0
                    matnam = 'AA';
                elseif lag == -1
                    matnam = 'BB';
                end
            end
            
            %--- yNames
            idx = indcv(vn, yNames);
            if idx > 0
                VarNr = idx;
                if lag == 0
                    matnam = 'CC';
                end
            end

            %--- zNames
            idx = indcv(vn, zNames);
            if idx > 0
                VarNr = idx;
                if lag == 0
                    matnam = 'DD';
                end
            end
            
        elseif strcmp(EqType, 'e');     %--- Expectational equations
            
            typ   = 'expectational';

            %--- xNames
            idx = indcv(vn, xNames);
            if idx > 0
                VarNr = idx;
                if lag == 1
                    matnam = 'FF';
                elseif lag == 0
                    matnam = 'GG';
                elseif lag == -1
                    matnam = 'HH';
                end
            end
            
            %--- yNames
            idx = indcv(vn, yNames);
            if idx > 0
                VarNr = idx;
                if lag == 1
                    matnam = 'JJ';
                elseif lag == 0
                    matnam = 'KK';
                end
            end

            %--- zNames
            idx = indcv(vn, zNames);
            if idx > 0
                VarNr = idx;
                if lag == 1
                    matnam = 'LL';
                elseif lag == 0
                    matnam = 'MM';
                end
            end
            
           
        else                            %--- Shock equation
            
            typ   = 'shock';

            %--- zNames
            idx = indcv(vn, zNames);
            if idx > 0
                VarNr = idx;
                if lag == 0
                    matnam = 'NN';  
                else
                    continue
                end
            else
                %--- Omit shock names
                continue
            end            
            
        end

        if strcmp(matnam, '');            
            error(['Variable ' varnames{i} ' not allowed in ' typ ' equation']);
        end
        
        linestr = [matnam '(' fillblanks(int2str(EqNr),2,'r') ',' fillblanks(int2str(VarNr),2,'r') ') = ' coeffs{i} '; '];
        linestr = fillblanks(linestr, SYSMAT_FIELDSIZE, 'l');
        
        str = [str linestr ' % ' varnames{i} crlf];
        
    end % for i = 1:length(varnames)

    
function idx = strfind2(str, pattern)
%**************************************************************************
% strfind2: Simulates behaviour of strfind using findstr 
%           to ensure compatibility with MATLAB version prior to 6.5.1
%**************************************************************************

if length(str) >= length(pattern)
    %--- String longer than search pattern --> perform search with findstr
    idx = findstr(str, pattern);
else
    %--- Search pattern longer than string to be searched -> return 0
    idx = [];
end


    
function idx = indcv(what,where)
%**************************************************************************
% Find character elements in cell array
%
% idx = indcv(what,where)
%
% what  [N,1] ... elements to be found
% where [M,K] ... cell array to be searched
%**************************************************************************

    if isempty(where)
        idx = 0;
        return
    end

    if isstr(what)
        what = cellstr(what);
    end
    
    if isstr(where)
        where = cellstr(where);
    end
    
    idx = zeros(size(what,1),size(where,2));

    for j=1:size(what,1)
        for c=1:size(where,2)
            for i=1:size(where,1)
                if strcmp(what{j},where{i,c}) == 1
                    idx(j,c) = i;
                    break;
                end
            end
        end
    end
      
    
function showdiagnostics(diagmat, VARIABLES, EQUATIONS, eqnames, eqtype, NrEqDet, NrEqExp, NrEqShk, NrVarState, NrVarJump, NrVarExo)
%**************************************************************************
% Show diagnostics for variables and equations and suggest types
%**************************************************************************

    %-- Count number of equations and variables
    NrEq     = NrEqDet     + NrEqExp   + NrEqShk;
    NrVarAll = NrVarState  + NrVarJump + NrVarExo;


    %--- Create list of shock names
    shklst = {};
    for i=1:size(VARIABLES,1)
        if strcmp(upper(VARIABLES(i,2)),'U') == 1
            shklst = [shklst; VARIABLES(i,1)];
        end
    end
    
    disp(' ')
    disp('========================================================================================================================================')
    disp('Your model consists of the following equations')
    disp('========================================================================================================================================')
    disp(' ')

    for i=1:size(eqnames,1)
        disp(eqnames{i})
        disp(EQUATIONS{i})
    end
    disp(' ')
    disp(' ')

    disp('========================================================================================================================================')
    disp('Diagnostics for equations')
    disp('========================================================================================================================================')
    disp(' ')   
    disp('For each equation, I will display the variables together with the time subscripts')
    disp('and will make some recommendations for the equation type.')
    disp('I use the following (heuristic) rules to suggest equation types (column ''Mine''):')   
    disp('- Equations with shock innovations (''U''; must be defined in VARIABLES) must be shock equations   --> ''S!''.')   
    disp('- Equations with leaded variables (t+1), which are not shock equations, must be expectational --> ''E!''.')   
    disp('- Other equations are suggested to be deterministic (''D''), but can also be expectational (''E'')--> ''D?''.')
    disp(' ')

    disp('------------------------------|Type -----|---------------------|-----------------------------------|-----------------------------------|');
    disp('Equation                      |Yours Mine|(t-1)                |(t)                                |(t+1)                              |');
    disp('------------------------------|---|------|---------------------|-----------------------------------|-----------------------------------|');

    
    varlag    = '';
    varcont   = '';
    varlead   = '';
    myeqtypes = cell(NrEq,1);

    for i=1:size(diagmat,1)
        
        varnam     = diagmat{i,1};
        lag        = diagmat{i,2};
        eqnr       = diagmat{i,3};
        if i < size(diagmat,1)
            eqnr_next = diagmat{i+1,3};
        else
            eqnr_next = 0;
        end

        if lag == 1
            varlead = [varlead varnam ', '];
        elseif lag  == 0
            varcont = [varcont varnam ', '];
        elseif lag  == -1
            varlag  = [varlag  varnam ', '];
        else
            err = 'You have an invalid lead/lag in equation: Stack the variable by ...';
        end
        
        if eqnr ~= eqnr_next
            %--- Next equation of end of list
            
              if indcv(varnam, shklst) 
                  mytype = 'S';
                  reliab = '!';
              elseif strcmp(varlead,'') == 0
                  mytype = 'E';
                  reliab = '!';
              else
                  mytype = 'D';
                  reliab = '?';
              end
              
            if length(varlead) > 2
                varlead = varlead(1:end-2);
            end
            if length(varcont) > 2
                varcont = varcont(1:end-2);
            end
            if length(varlag) > 2
                varlag = varlag(1:end-2);
            end

            myeqtypes{eqnr} = mytype;
            
            eqnam    = eqnames{eqnr, 1};
            eqnam    = eqnam(1:min(30,length(eqnam)));
            yourtype = eqtype{eqnr, 1};
            disp([fillblanks(eqnam,30,'l') '| ' yourtype ' | ' fillblanks([mytype reliab],5,'l') '| '  ... 
                  fillblanks(varlag,20,'l') '|' fillblanks([fillblanks(varcont,35,'l') '|' varlead],71,'l') '|'])

            varlag  = '';
            varcont = '';
            varlead = '';    
        end
    end
    disp('------------------------------|---|------|---------------------|-----------------------------------|-----------------------------------|');
    disp(' ')   
    disp(' ')   
    
  



    disp('========================================================================================================================================')
    disp('Diagnostics for variables (which occure in equations)')
    disp('========================================================================================================================================')
    disp(' ')
    disp('For each variable, I will display the equation numbers in which they occur, together with the time subscript')
    disp('and will make some recommendations for the variable type.')
    disp('I use the following (heuristic) rules to suggest variable types (column ''Mine''):')   
    disp('- You have selected the variable to be a shock innovation (''U''): I trust you --> ''U!''')   
    disp('- Variable occurs in shock equation (''S''): Variable must be an exogenous state variable --> ''Z!''')   
    disp('- Variable occurs at t-1, t and t+1: Variable must be an endogenous state variable --> ''X!''')   
    disp('- Variable occurs at t-1: Variable must be an endogenous state variable (but can in principle be shifted in time) --> ''X''')   
    disp('- Other variables can either be endogenous state variables (''X'') or other endogenous (jump) variables (''Y'') --> ''X/Y?''')   
    disp(' ')   
    disp(' ')   
    
    diagmatvar = sortrows(diagmat, 1);
    varnam_next = '';
    eqlag     = '';
    eqcont    = '';
    eqlead    = '';
    eqerr     = '';
    eqerrtxt  = '';
    

    disp('------------------------------|Type -----|---------------------|-----------------------------------|-----------------------------------|');
    disp('Variable                      |Yours Mine|(t-1)                |(t)                                |(t+1)                              |');
    disp('------------------------------|---|------|---------------------|-----------------------------------|-----------------------------------|');

    
    for i=1:size(diagmatvar,1)
        
        varnam     = diagmatvar{i,1};
        lag        = diagmatvar{i,2};
        eqnr       = diagmatvar{i,3};
        op         = diagmatvar{i,4};
        
        if i < size(diagmatvar,1)
            varnam_next = diagmatvar{i+1,1};
        else
            varnam_next = '';
        end

        if lag == 1
            eqlead = [eqlead int2str(eqnr) ', '];
        elseif lag  == 0
            eqcont = [eqcont int2str(eqnr) ', '];
        elseif lag  == -1
            eqlag  = [eqlag  int2str(eqnr) ', '];
        else
            eqerrtxt = [eqerrtxt 'Variable ''' op ''' has invalid lag. Stack the variable by including an addition equation..' char(13) ]
        end
        
        if strcmp(varnam, varnam_next) == 0
            %--- Next variable of end of list

            idx = indcv(varnam,VARIABLES(:,1));
            if idx > 0
                yourtype = VARIABLES{idx,2};
            else
                yourtype = '-';
            end

            if strcmp(eqlead,'') == 0 & strcmp(eqcont,'') == 0 & strcmp(eqlag,'') == 0
                %--- (t-1), (t) and (t+1) --> 'X'
                mytype = 'X!  ';
            elseif strcmp(eqlag,'') == 0
                %--- (t-1) --> 'X'
                mytype = 'X   ';
            elseif strcmp(yourtype, 'U')
                mytype = 'U!  ';
            elseif strcmp(myeqtypes(eqnr),'S')
                %--- If variable occurs in shock equation ('S') --> 'Z'
                mytype = 'Z!  ';
            else
                mytype = 'X/Y?';
            end
            
            if length(eqlead) > 2
                eqlead = eqlead(1:end-2);
            end
            if length(eqcont) > 2
                eqcont = eqcont(1:end-2);
            end
            if length(eqlag) > 2
                eqlag = eqlag(1:end-2);
            end
            
            disp([fillblanks(varnam,30,'l') '| ' yourtype ' | ' mytype ' | ' fillblanks(eqlag,20,'l') '|' fillblanks(eqcont,35,'l') '|' fillblanks(eqlead,35,'l') '|'])

            
            eqlag  = '';
            eqcont = '';
            eqlead = '';
        end
           
    end

    disp('------------------------------|---|------|---------------------|-----------------------------------|-----------------------------------|');
    disp(' ')   
    disp(eqerrtxt )
    disp(' ')   

    
    if NrEq == NrVarAll 
        disp(['Your model consists of ' int2str(NrEq) ' equations and ' int2str(NrVarAll) ' variables. This is fine!'])
    elseif NrEq > NrVarAll
        error(['Your model consists of ' int2str(NrEq) ' equations, but only of ' int2str(NrVarAll) ' variables. This won''t make you happy !!!!!!!!'])
    else
        error(['Your model consists of ' int2str(NrVarAll) ' variables (as defined in ''VARIABLES''), but only of ' int2str(NrEq) ' equations. This won''t make you happy !!!!!!!!'])
    end
    disp(['   thereof ' int2str(NrEqDet) ' deterministic equation(s),'])
    disp(['           ' int2str(NrEqExp) ' expectational equation(s),'])
    disp(['           ' int2str(NrEqShk) ' shock equation(s),'])
    disp(' ')   
    disp(['           ' int2str(NrVarState) ' endogenous state variable(s),'])
    disp(['           ' int2str(NrVarJump) ' endogenous jump variable(s) and'])
    disp(['           ' int2str(NrVarExo) ' exogenous state variable(s).'])
    disp('    Remark: I do not count innovations to shocks (''U'').')   
    disp(' ')   
    disp(' ')   
    
    if NrVarJump > NrEqDet
        error('You have more jump variables than deterministic equations. This causes problems with the C matrix !!!!!!!!')   
    else
        disp('Your number of deterministic equations is >= than your number of jump variables.')
        disp('This is fine, since the dimensions of matrix CC are OK.') 
        disp('Neverteless, please checke whether CC has full rank (especially that there are no zero rows or columns.')
        disp('If this is not the case, redeclare (one ore more) jump variables to be state variable(s).')   
    end
    disp(' ')   
    disp('If you are not sure about the variable and equation types you can solve your model with brute force by declaring all')   
    disp('endogenous variables as state variables (''X'') and all equations (w.o. the shock equations) as expectational equations (''E'').')
    disp('This is inefficient and may cause troubles with larger models, but it might help you to get your model running. ')   
    disp(' ') 
    disp('You can turn off the diagnostics by setting SYSMAT_DISPDIAG=0')
    disp(' ')   
    disp(' ') 
    