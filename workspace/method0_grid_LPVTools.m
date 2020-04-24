function [ret] = method0_grid_LPVTools(modelname, A_fh, B_fh, C_fh, D_fh, p_lim, dp_lim, bases, bases_Jac, Resolution)
%% LPV_L2anal_LPVTools
%
%  File: LPV_L2anal_LPVTools.m
%  Directory: 1_PhD_projects/22_Hinf_norm/LPV_TAC_v2_newmod (original dir)
%  Author: Peter Polcz (ppolcz@gmail.com)
%
%  Created on 2018. November 19.
%  Major review on 2020. March 17. (2019b)
%

%%

TMP_uBUOdotSvmFyUWQKMUdr = pcz_dispFunctionName('LPVTools');

%% Gridded model

np = size(p_lim,1);

p = sym('p',[np 1]);

p_cell = cell(1,np);
p_name = cellfun(@(i) {['p' num2str(i)]}, num2cell(1:np));
for i = 1:np
    p_cell{i} = pgrid(p_name{i}, linspace(p_lim(i,1),p_lim(i,2),Resolution(i)), dp_lim(i,:));
end

[A_unc,B_unc,C_unc,D_unc] = helper_convert(A_fh,B_fh,C_fh,D_fh,p_cell,'pmat');

% Convert symbolic basis functions to pmat
symobj = num2cell([bases,bases_Jac]);
fhobj = cellfun(@(obj) { matlabFunction(obj,'vars',p) }, symobj );
pmatobj = cellfun(@(obj) { obj(p_cell{:}) }, fhobj);

% Define basis functions:
args = cellfun(@(fi,gradi) {[ {fi} reshape([p_name ; gradi],[1 2*np])]}, pmatobj(:,1), num2cell(pmatobj(:,2:end),2));
bases = cellfun(@(args) basis(args{:}), args);

sys_unc_GRID = ss(A_unc,B_unc,C_unc,D_unc);

%% Computations

TMP_cwCXkgHfZmFQRzNVUlCO = pcz_dispFunctionName('LPVTools grid based', 'lpvnorm');

[ lpv_gamma, X, info ] = lpvnorm(sys_unc_GRID, bases);
Overall_Time = toc(TMP_cwCXkgHfZmFQRzNVUlCO);

pcz_dispFunction('Solver time: <strong>%g</strong>', Overall_Time)
pcz_dispFunction(2, sprintf('<strong>gamma = %g </strong>(lpvnorm)', lpv_gamma))

pcz_dispFunctionEnd(TMP_cwCXkgHfZmFQRzNVUlCO);
% -------------------------------------------------------------------------

info = sprintf('pgrid (%s)x(%s)', ...
    strjoin(cellfun(@(n) {num2str(n)}, num2cell(Resolution) ), 'x'),...
    strjoin(cellfun(@(n) {num2str(n)}, num2cell(0*Resolution + 2) ), 'x'));

store_results('Results_All', modelname, 0, lpv_gamma, 0, Overall_Time, info, 'LPVTools - grid')

%%

pcz_dispFunctionEnd(TMP_uBUOdotSvmFyUWQKMUdr);

end
