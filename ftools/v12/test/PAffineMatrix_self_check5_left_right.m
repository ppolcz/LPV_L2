%% PAffineMatrix_self_check1_left_right
%  
%  File: PAffineMatrix_self_check1_left_right.m
%  Directory: 7_ftools/ftools/v11/test
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2019. November 05. (2019a)
%

syms p real

xi = [1;p];

N = PAffineMatrix([1 2 3 4 5 6], xi, p, 'type', 'left')

Nr_sym = sym(N)

Theta_right = N.Theta

N.type = PAffineMatrix.TYPE_LEFT

Theta_left = N.Theta

s = N.s

m = N.m

Nl_sym = sym(N)



%  ITT HAGYTAD ABBA:
% 
% Nl_sym    NEM EGYENLO     Nr_sym-el

