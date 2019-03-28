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
    
