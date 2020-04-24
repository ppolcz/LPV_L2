%% PAffineMatrix_self_check6_adapt_channels
%  
%  File: PAffineMatrix_self_check6_adapt_channels.m
%  Directory: 7_ftools/ftools/v11/test
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2019. November 05. (2019a)
%

%%

TMP_LiCDAYlkjszFUeQxehlm = pcz_dispFunctionName;

syms x p r real

X = PAffineMatrix([1 2 3 4 5 6],[1 x],'name','X');

Y = PAffineMatrix([1 2 3 4 5 6 7 8 9],[p r 1],'name','Y');

% 2019.11.05. (november  5, kedd), 20:16 [NEM MUKODIK]
[X,Y] = adapt_channels(X,Y)


pcz_dispFunctionEnd(TMP_LiCDAYlkjszFUeQxehlm);