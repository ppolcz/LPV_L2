function [ret] = method3_IQC_LFT_LPVTools(modelname,A_fh,B_fh,C_fh,D_fh,p_lim,dp_lim)
%% LPV_L2anal_LPVTools
%
%  File: LPV_L2anal_LPVTools.m
%  Directory: 1_PhD_projects/22_Hinf_norm/LPV_TAC_v2_newmod (original dir)
%  Author: Peter Polcz (ppolcz@gmail.com)
%
%  Created on 2018. November 19.
%  Major review on 2020. March 17. (2019b)
%
% ``IQCs are used for worst-case analysis (lpvwcgain) of uncertain LPV
% systems (grid-based [3] and LFT-based [4,5,6]), and for analysis
% (lpvnorm) of nominal (not uncertain) rate-bounded LFT-based LPV systems
% [7].'' (LPVToolbox Manual)
% 
% [3] Pfifer and Seiler (2015). Robustness analysis of linear parameter
% varying systems using integral quadratic constraints. International
% Journal of Robust and Nonlinear Control, Wiley Online Library,
% 25(15):2843-2864, 2015.
% 
% [4] C. Scherer and S. Wieland, "Linear matrix inequalities in control,"
% Lecture notes for a course of the dutch institute of systems and control,
% Delft University of Technology, 2004.
% 
% [5] C. Scherer and I. Kose, "Robustness with dynamic IQCs: An exact
% state-space characterization of nominal stability with applications to
% robust estimation," Automatica, Vol. 44, No. 7, pp. 1666-1675, 2008.
% 
% [6] C. Scherer, "LPV control and full-block multipliers," Automatica,
% Vol. 37, No. 3, pp. 361-375, 2001.
% 
% [7] Helmersson (1999). An IQC-based stability criterion for systems with
% slowly varying parameters. IFAC Proceedings Volumes, 32(2):3183-3188,
% 1999.


%%

TMP_uBUOdotSvmFyUWQKMUdr = pcz_dispFunctionName('LPVTools');

%% LFT model

np = size(p_lim,1);
p_cell = cell(1,np);
for i = 1:np
    name = [ 'p' num2str(i) ];
    p_cell{i} = tvreal(name, p_lim(i,:), dp_lim(i,:));
end

[A_unc,B_unc,C_unc,D_unc] = helper_convert(A_fh,B_fh,C_fh,D_fh,p_cell,'plftmat');

sys_unc_LFT = ss(A_unc,B_unc,C_unc,D_unc);
simplify(sys_unc_LFT,'full');

%% Computations -- LPVNORM

% -------------------------------------------------------------------------
TMP_cwCXkgHfZmFQRzNVUlCO = pcz_dispFunctionName('LPVTools LFT/IQC', 'lpvnorm');

% What actually is executed:
% 
% [7] Helmersson (1999). An IQC-based stability criterion for systems with
% slowly varying parameters. IFAC Proceedings Volumes, 32(2):3183-3188,
% 1999.
lpv_gamma = lpvnorm(sys_unc_LFT);
lpvnorm_time = toc(TMP_cwCXkgHfZmFQRzNVUlCO);

pcz_dispFunction('Solver time: <strong>%g</strong>', lpvnorm_time)
pcz_dispFunction(2, sprintf('<strong>gamma = %g </strong>(lpvnorm)', lpv_gamma))

pcz_dispFunctionEnd(TMP_cwCXkgHfZmFQRzNVUlCO);
% -------------------------------------------------------------------------

store_results('LPVTools_Results.csv', modelname, 0, lpv_gamma, 0, lpvnorm_time, 'tvreal, lpvnorm', 'LFT/IQC')
store_results('Results_All.csv', modelname, 0, lpv_gamma, 0, lpvnorm_time, 'tvreal, lpvnorm', 'LPVTools - LFT/IQC')

%% Computations -- LPVWCGAIN

% -------------------------------------------------------------------------
TMP_KZWeXiYFmdpQdsgidKeG = pcz_dispFunctionName('LPVTools LFT/IQC', 'lpvwcgain');

% What actually is executed:
% 
% [4] C. Scherer and S. Wieland, "Linear matrix inequalities in control,"
% Lecture notes for a course of the dutch institute of systems and control,
% Delft University of Technology, 2004.
% 
% [5] C. Scherer and I. Kose, "Robustness with dynamic IQCs: An exact
% state-space characterization of nominal stability with applications to
% robust estimation," Automatica, Vol. 44, No. 7, pp. 1666-1675, 2008.
% 
% [6] C. Scherer, "LPV control and full-block multipliers," Automatica,
% Vol. 37, No. 3, pp. 361-375, 2001.
wcg = lpvwcgain(sys_unc_LFT);
lpvwcgain_time = toc(TMP_KZWeXiYFmdpQdsgidKeG);

bounds = zeros(1,2);
if wcg.LowerBound
    bounds(1) = wcg.LowerBound;
end
if wcg.UpperBound
    bounds(2) = wcg.UpperBound;
end

pcz_dispFunction('Solver time: <strong>%g</strong>', lpvwcgain_time)
pcz_dispFunction(2, 'Bounds: [<strong>%g</strong>,%g] = [%g,%g] dB ', bounds, 20*log10(bounds))

pcz_dispFunctionEnd(TMP_KZWeXiYFmdpQdsgidKeG);
% -------------------------------------------------------------------------

store_results('LPVTools_Results.csv', modelname, bounds(1), bounds(2), 0, lpvwcgain_time, 'tvreal, lpvwcgain', 'LFT/IQC')
store_results('Results_All.csv', modelname, bounds(1), bounds(2), 0, lpvwcgain_time, 'tvreal, lpvwcgain', 'LPVTools - LFT/IQC')

%%

pcz_dispFunctionEnd(TMP_uBUOdotSvmFyUWQKMUdr);

end
