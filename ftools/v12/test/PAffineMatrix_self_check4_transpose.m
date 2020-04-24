%% 
%  
%  File: PAffineMatrix_self_check4_transpose.m
%  Directory: 1_PhD_projects/00_my_toolboxes/FinslerTools/v11/test
%  Author: Peter Polcz (ppolcz@gmail.com) 
% 
%  Created on 2018. October 23.

%%

syms x1 x2 real
x = [x1;x2];

Theta = [1 0 0 1 ; -1 1 1 0];

N = PAffineMatrix(Theta,x);
N.Sym = 1;

Nct = N';
Nt = N.';

pcz_symzero(N.Sym' - Nct.Sym,'N(x) conjucate transpose')
pcz_symzero(N.Sym.' - Nt.Sym,'N(x) transpose')
