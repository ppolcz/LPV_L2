function [N] = to_symTheta(N,name,assumptions)
%% to_symTheta
%  
%  File: to_symTheta.m
%  Directory: 1_PhD_projects/00_my_toolboxes/algo_P/@PGenAffineMatrix
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  TODO
%

%%

if nargin > 1
    N.name = name;
end

if nargin <= 2
    assumptions = 'real';
end

if isa(N.Theta,'sdpvar')
    Theta = pcz_sdpvar2sym(N.name,N.Theta,assumptions);

    % svar = sym(N.Theta);
    % variables = symvar(svar);
    % new_vars = sym(name, size(variables), assumption);
    % ret = subs(svar, variables, new_vars);

else
    Theta = sym(N.Theta);
end

if isa(N.channels,'sdpvar')
    [channels,vars,new_vars] = pcz_sdpvar2sym(N.name,N.channels,assumptions);
    subsvars = subs(N.subsvars,vars,new_vars);
else
    channels = sym(N.channels);
    subsvars = sym(N.subsvars);
end

N = PGenAffineMatrix(Theta, channels, subsvars);
N = N.generate_symbolic;

end