function [varargout] = helper_fh2fh_vec(A_fh, B_fh, C_fh, D_fh, p_lim)
%% helper_fh2fh_vec
%
%  File: helper_fh2fh_vec.m
%  Directory: 8_published/LPV_L2/workspace
%  Author: Peter Polcz (ppolcz@gmail.com)
%
%  Created on 2020. March 25. (2019b)
%

%%

p = sym('p', [size(p_lim,1),1]);

symobj = cell(1,7);

[symobj{:}] = helper_fh2sym(A_fh,B_fh,C_fh,D_fh,p_lim);

varargout = cellfun( @(obj) { matlabFunction(obj,'vars', {p}) }, symobj );

end