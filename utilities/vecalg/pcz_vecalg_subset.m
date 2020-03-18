function [ret] = pcz_vecalg_subset(V, W, tol)
%% Script pcz_vecalg_subset
%  
%  File: pcz_vecalg_subset.m
%  Directory: projects/3_outsel/2018_01_10_LPV_inversion
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2018. January 13.
%
%%

if nargin < 3
    tol = 1e-10;
end

W = pcz_vecalg_orth(W,tol);
ret = rank(W / (W'*W) * W' * V - V, tol) == 0;

end