function [Ker_X] = pcz_vecalg_null(X, tol)
%% Script pcz_vecal_null
%  
%  File: pcz_vecalg_null.m
%  Directory: projects/3_outsel/2018_01_10_LPV_inversion
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2018. January 13.
%
% tol: rank tolerance
% X: the same as Matlab's build in |null|
% X = [ v_1
%       v_2
%       ...
%       v_n ], where v_i are row vectors.
% 
% Return: Ker_X, such that X * Ker_X = 0
% 

%%
% To be as implemented in Matlab's build-in |null|
X = X';

if nargin < 2
    tol = 1e-10;
end

%%

Im_X = pcz_vecalg_orth(X,tol);
Ker_X = null(Im_X');

assert(norm(Ker_X' * X) < tol, ...
    'Ker_X and Im_X are not orthogonal. Rank tolerance does not fulfill.');

assert(rank([X , Ker_X],tol) == size(X,1), ...
    'Ker_X + Im_X are not complete. Rank tolerance does not fulfill.')

end