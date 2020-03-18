function [ret,variables,new_vars] = pcz_sdpvar2sym(name, var, assumption)
%% pcz_sdpvar2sym
%  
%  File: pcz_sdpvar2sym.m
%  Directory: 2_demonstrations/lib/matlab
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2018. September 17.
%
% Example (given in 2019.11.05. (november  5, kedd), 20:25):
% 
%  >> [ret,var,newvar] = pcz_sdpvar2sym('P',sdpvar(2),'real')
%   
%  ret = [ P1, P2]
%        [ P2, P3]
%  
%  var = [ x4, x5, x6]
%  
%  newvar = [ P1, P2, P3]

%%

if nargin <= 2
    assumption = 'real';
end

svar = sym(var);

variables = symvar(svar);

new_vars = sym(name, size(variables), assumption);

ret = subs(svar, variables, new_vars);

end