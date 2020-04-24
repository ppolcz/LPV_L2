function [A_sym,B_sym,C_sym,D_sym,AC_sym,BD_sym,M_sym] = helper_fh2sym(A_fh, B_fh, C_fh, D_fh, p_lim)
%% helper_fh2sym
%
%  File: helper_fh2sym.m
%  Directory: 8_published/LPV_L2/workspace
%  Author: Peter Polcz (ppolcz@gmail.com)
%
%  Created on 2020. March 25. (2019b)
%

%%

p = sym('p',[size(p_lim,1),1]);
p_cell = num2cell(p);

% Symbolical expression of matrix A(p):
if isnumeric(A_fh)
    A_sym = sym(A_fh);
else
    A_sym = A_fh(p_cell{:});
end

% Symbolical expression of matrix B(p):
if isnumeric(B_fh)
    B_sym = sym(B_fh);
else
    B_sym = B_fh(p_cell{:});
end

% Symbolical expression of matrix C(p):
if isnumeric(C_fh)
    C_sym = sym(C_fh);
else
    C_sym = C_fh(p_cell{:});
end

% Symbolical expression of matrix D(p):
if isnumeric(D_fh)
    D_sym = sym(D_fh);
else
    D_sym = D_fh(p_cell{:});
end

AC_sym = [
    A_sym
    C_sym
    ];

BD_sym = [
    B_sym
    D_sym
    ];

M_sym = [ AC_sym BD_sym ];

end