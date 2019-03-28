%--- Written by Martin Schneider, Oesterreichische Nationalbank, 2007 (martin.schneider@oenb.at)
%--- Version 1.0: July 2007

The MATLAB program makesysmat greatly simplifies the way how to write down and solve a model with Uhlig's Toolkit. You just have to write down the loglinear equations of your model in a convenient form (e.g. '0 = - y(t) + z(t) + rho * k(t-1)') and specify the type of each equation and variable. Then, makesysmat does the cumbersome work of translating the model into the matrices AA-NN for you. It carefully parses the equations for syntax errors and provides a diagnostic tool that helps you to properly set the equation and variable types. The program writes MATLAB code into the file "filename_SYSMAT.m"and runs that code in one step. 

You can easily integrate the program into your models. You just have to copy the files callsysmat.m and makesysmat.m in the Toolkit folder and add "callsysmat" as the first command in the Toolkit file do_it.m. The program automatically recognizes whether you have a usual Toolkit m-file (with matrices AA-NN specified) or an m-file that uses makesysmat.m. 

Run exampl0s.m to see how the program works and try to create your own version of the model!
