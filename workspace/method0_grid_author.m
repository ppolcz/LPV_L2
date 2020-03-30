function method0_grid_author(modelname, A_fh, B_fh, C_fh, D_fh, p_lim, dp_lim, bases, bases_Jac, Resolution)
%%
%  File: method0_grid_author.m
%  Directory: 8_published/LPV_L2/workspace
%  Author: Tamas Peni, Peter Polcz (ppolcz@gmail.com) 
% 
%  Major review on 2020. March 26. (2019b)

[A_fh,B_fh,C_fh,D_fh] = helper_fh2fh_vec(A_fh,B_fh,C_fh,D_fh,p_lim);

pargrd = [ cellfun(@(c) {linspace(c{:})}, num2cell(num2cell([p_lim Resolution']),2))' num2cell(dp_lim,2)'];

np = size(p_lim,1);

p = sym('p',[np 1]);
dp = sym('dp',[np 1]);
pdp = [ p ; dp ];

Pbase = matlabFunction(bases, 'vars', {pdp.'});
Pbase_der = matlabFunction(bases_Jac*dp, 'vars', {pdp.'});

OVERALL = tic;
[gamma,Pvars_grid_Tamas,Solver_Time] = lpvL2gain(@lpvsys,pargrd,Pbase,Pbase_der);

pcz_2basews(Pvars_grid_Tamas)

info = sprintf('(%s)x(%s)', ...
    strjoin(cellfun(@(n) {num2str(n)}, num2cell(Resolution) ), 'x'),...
    strjoin(cellfun(@(n) {num2str(n)}, num2cell(0*Resolution + 2) ), 'x'));

store_results('Results_All.csv', modelname, 0, gamma, Solver_Time, toc(OVERALL), ...
    info, 'Tamas - grid')

function [Ap,Bp,Cp,Dp,S]=lpvsys(p)
    Ap = A_fh(p');
    Bp = B_fh(p');
    Cp = C_fh(p');
    Dp = D_fh(p');
end

end