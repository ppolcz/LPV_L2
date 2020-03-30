function method1_grid_uncertain_LPVTools(modelname, A_fh, B_fh, C_fh, D_fh, p_lim, dp_lim, Resolution)
%% 
%  File: method3_grid_uncertain_LPVTools.m
%  Directory: 8_published/LPV_L2/workspace
%  Author: Peter Polcz (ppolcz@gmail.com) 
% 
%  Created on 2020. March 26. (2019b)

% ``IQCs are used for worst-case analysis (lpvwcgain) of uncertain LPV
% systems (grid-based [3] and LFT-based [*]), (...)'' (LPVToolbox Manual)
% 
% [3] Pfifer and Seiler (2015). Robustness analysis of linear parameter
% varying systems using integral quadratic constraints. International
% Journal of Robust and Nonlinear Control, Wiley Online Library,
% 25(15):2843-2864, 2015.

%%

TMP_uBUOdotSvmFyUWQKMUdr = pcz_dispFunctionName('LPVTools');

%% Gridded model

np = size(p_lim,1);
p_nom = sum(p_lim,2)/2;

p_cell = cell(1,np);
p_cell{1} = pgrid('p1',linspace(p_lim(1,1),p_lim(1,2),Resolution(1)), dp_lim(1,:));
for i = 2:np
    name = [ 'p' num2str(i) ];
    p_cell{i} = ureal(name, p_nom(i), 'Range', p_lim(i,:));
end

A_fh2 = @(p1,p2,p3) [
    35/(3*(p2^2 - 7)), (10*p2)/(p2^2 - 7), (7*p3)/(3*(p2^2 - 7))
                    0,                  0,                     1
    (5*p2)/(p2^2 - 7),      30/(p2^2 - 7),    (p2*p3)/(p2^2 - 7)
    ];

[A_unc2,B_unc,C_unc,D_unc] = helper_fh2lfr(A_fh2,B_fh,C_fh,D_fh,p_cell,'umat');

A_unc = A_unc2 * blkdiag(1,p_cell{1},1);

sys_unc = ss(A_unc,B_unc,C_unc,D_unc);

%% Computations

TMP_cwCXkgHfZmFQRzNVUlCO = pcz_dispFunctionName('LPVTools LTI worst case analysis', 'lpvwcgain');

[ lpv_gamma, X, info ] = lpvwcgain(sys_unc);
Overall_Time = toc(TMP_cwCXkgHfZmFQRzNVUlCO);

pcz_dispFunction('Solver time: <strong>%g</strong>', Overall_Time)
pcz_dispFunction(2, sprintf('<strong>gamma = %g </strong>(lpvnorm)', lpv_gamma))

pcz_dispFunctionEnd(TMP_cwCXkgHfZmFQRzNVUlCO);
% -------------------------------------------------------------------------

info = sprintf('pgrid (%s)x(%s)', ...
    strjoin(cellfun(@(n) {num2str(n)}, num2cell(Resolution) ), 'x'),...
    strjoin(cellfun(@(n) {num2str(n)}, num2cell(0*Resolution + 2) ), 'x'));

store_results('LPVTools_Results.csv', modelname, 0, lpv_gamma, 0, Overall_Time, info, 'grid')
store_results('Results_All.csv', modelname, 0, lpv_gamma, 0, Overall_Time, info, 'LPVTools - grid')

%%

pcz_dispFunctionEnd(TMP_uBUOdotSvmFyUWQKMUdr);

end
