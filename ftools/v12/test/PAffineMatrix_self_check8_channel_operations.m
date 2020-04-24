%% PAffineMatrix_self_check8_channel_operations
%  
%  File: PAffineMatrix_self_check8_channel_operations.m
%  Directory: 7_ftools/ftools/v11/test
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2019. November 05. (2019a)
%

%% set_channel_by_index

syms a b c real
N = PAffineMatrix(zeros(4,12),[1;a;b;c],'type','right','name','N');
N = N.set_channel_by_index(ones(4,3),2)
N_sym = N.Sym

M = PAffineMatrix(zeros(4,12),[1;a;b;c],'type','left','name','M');
M = M.set_channel_by_index(ones(4,3),2)
M_Sym = M.Sym

M_plfr = plfr(M);

pcz_symzero(sym(M_plfr) - M_Sym, 'plfr = PAffineMatrix')

%%


