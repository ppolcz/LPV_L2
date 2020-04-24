function [S, Pi_hat, S_plus] = P_get_basis_rationals(Pi, w)
%% P_get_basis_rationals
%  
%  File: P_get_basis_rationals.m
%  Directory: 1_PhD_projects/00_my_toolboxes/FinslerTools/v11
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2018. November 19.
%

%% Initialization

[Theta,~,~] = P_Pi_canonical_decomp(Pi,w);

[m,K] = size(Theta);
mh = rank(Theta);
% pcz_display(n,p,k,K);

%% Permutation of variables of LFR
% Permutation matrix ${\tt Ir} = I_\varrho$.
[~,Icols] = rref(Theta);
Jcols = setdiff(1:K,Icols);
rho = [Icols Jcols];
Ir = pcz_permat(rho)';

%%
% Permutation matrix ${\tt Isb} = I_{\sigma_b}$ and ${\tt Is} = I_\sigma$.
[~,Irows] = rref(Theta');
Jrows = setdiff(1:m,Irows);
sigma_b = [Irows Jrows];
Isb = pcz_permat(sigma_b)';
% Is = Isb(n+1:end,n+1:end);

% pcz_display(Irows,Jrows,sigma_b,Isb,Is)

%%
% Permuted coefficient matrix
pTheta = Theta(sigma_b,rho);

% pcz_display(pTheta)

assert(rank(pTheta(1:mh,1:mh)) == mh, ...
    'The upper left (mh)x(mh) submatrix of Theta should be full-rank');

%% 
% Auxiliary matrices
V = pTheta(1:mh,:);
W = pTheta(mh+1:end,:);
Gamma = round(W*V'/(V*V'),10);

S = [ eye(mh) ; Gamma ];
S(sigma_b,:) = S;

S_plus = pinv(S);

Pi_hat = S_plus * Pi;


end
