function [X] = to_symTheta(X,name,assumptions)
%% to_symTheta
%  
%  File: to_symTheta.m
%  Directory: 1_PhD_projects/00_my_toolboxes/algo_P/@PAffineMatrix
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2018. October 21.
%

%%

if nargin > 1
    X.name = name;
end

if nargin <= 2
    assumptions = 'real';
end

if isa(X.Theta,'sdpvar')
    Theta = pcz_sdpvar2sym(X.name,X.Theta,assumptions);
else
    Theta = sym(X.Theta);
end

if isa(X.channels,'sdpvar')
    [channels,vars,new_vars] = pcz_sdpvar2sym(X.name,X.channels,assumptions);
    subsvars = subs(X.subsvars,vars,new_vars);
else
    channels = sym(X.channels);
    subsvars = sym(X.subsvars);
end

X.channels = channels;
X.subsvars = subsvars;
X.Theta = Theta;
X.Sym = 1;

% X = PAffineMatrix(Theta, channels, subsvars);
% X = copy(X,X);
% X.Sym = 1;

end