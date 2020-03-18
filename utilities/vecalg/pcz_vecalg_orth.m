function [Im_X] = pcz_vecalg_orth(X, tol)
%% Script pcz_vecalg_orth
%  
%  File: pcz_vecalg_orth.m
%  Directory: demonstrations/lib/matlab
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2018. January 13.
%
%%

if nargin < 2
    tol = 1e-10;
end

%%

[U,S,~] = svd(X);

Im_X = U(:,1:rank(S,tol));


end