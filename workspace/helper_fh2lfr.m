function [A_lfr,B_lfr,C_lfr,D_lfr,AC_lfr,BD_lfr,M_lfr] = helper_fh2lfr(A_fh, B_fh, C_fh, D_fh, p_lim)
%% helper_fh2lfr
%  
%  File: helper_fh2lfr.m
%  Directory: 8_published/LPV_L2/workspace
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2020. March 25. (2019b)
%

%%

[~,p_cell] = pcz_generateLFRStateVector('p',p_lim);

% LFR realization of matrix A(p):
if isdouble(A_fh)
    A_lfr = lfr(A_fh);
else
    A_lfr = A_fh(p_cell{:});
end

% LFR realization of matrix B(p):
if isdouble(B_fh)
    B_lfr = lfr(B_fh);
else
    B_lfr = B_fh(p_cell{:});
end

% LFR realization of matrix C(p):
if isdouble(C_fh)
    C_lfr = lfr(C_fh);
else
    C_lfr = C_fh(p_cell{:});
end

% LFR realization of matrix D(p):
if isdouble(D_fh)
    D_lfr = lfr(D_fh);
else
    D_lfr = D_fh(p_cell{:});
end

AC_lfr = [
    A_lfr
    C_lfr
    ];

BD_lfr = [
    B_lfr
    D_lfr
    ];

M_lfr = [ AC_lfr BD_lfr ];

end