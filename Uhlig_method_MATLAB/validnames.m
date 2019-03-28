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
    