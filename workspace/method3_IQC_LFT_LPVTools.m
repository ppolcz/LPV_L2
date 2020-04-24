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
% ``In the LFT case, lpvnorm implements a computation of the induced L2
% norm only based on the derivation by [1]**. lpvwcgain computes the upper
% bound on the induced L2 norm of an uncertain LPV system (grid- or
% LFT-based). Its implementation is based on [3], and [8]***.'' ([2])
% 
% 
% ** Author's comment: when the rate bounds are not specified. Otherwise,
% `plftss/lpvnorm' calls `plftss/lpvwcgain'.
% 
% *** Author's comment: The structure of the dynamic multiplier block
% corresponding to the smooth time-varying parametric uncertainty is
% similar to that in [7] obtained through the swapping lemma. Differently
% from [7], in which only an LPV stability test is proposed,
% lpvnorm/lpvwcgain can handle a performance channel (multiplier block) to
% compute an upper-bound on the L2-gain.
% 
% 
% 1. In lpvwcgain, the basis functions for the dynamic multiplier are
% considered similarly to [8,9]. The difference is that instead of the
% repeated poles [ 1 1/(s-p) ... 1/(s-p)^k ] of [8,9], lpvnorm/lpvwcgain
% considered multiple poles [1 1/(s-p1) ... 1/(s-pk) ] appearing once in a
% multiplier block.
% 
% 2. lpvwcgain does NOT solve the KYP LMI associated with the IQC analysis.
% Instead, it enforces the corresponding frequency domain constraint on a
% frequency grid. The implementation is iterative based on the cutting
% plane algorithm of [10].
% 
% 
% References
% 
% [1] Apkarian and Gahinet (1995). A convex characterization of
% gain-scheduled H_∞ controllers. IEEE Transactions on Automatic Control,
% 40(5):853-864, 1995.
% 
% [2] Hjartarson, Seiler and Packard (2015). LPVTools: a toolbox for
% modeling, analysis, and synthesis of parameter varying control systems.
% IFAC-PapersOnLine, 48(26):139-145, 2015.
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
% 
% [8] Veenman and Scherer (2014). IQC-synthesis with general dynamic
% multipliers. International Journal of Robust and Nonlinear Control, Wiley
% Online Library, 24(17):3027-3056, 2014.
% 
% [9] Köroğlu and Scherer (2006). Robust stability analysis against
% perturbations of smoothly time-varying parameters. 45th IEEE Conference
% on Decision and Control, :2895-2900, 2006.
% 
% [10] Wallin, Kao and Hansson (2008). A cutting plane method for solving
% KYP-SDPs. Automatica, 44(2):418-429, 2008.
% 

global LPVTools_vars
LPVTools_vars.Solver_Time = 0;

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
sys_unc_LFT = simplify(sys_unc_LFT,'full');

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

store_results('Results_All', modelname, bounds(1), bounds(2), LPVTools_vars.Solver_Time, lpvwcgain_time, 'tvreal, lpvwcgain', 'LPVTools - LFT/IQC')

%%

pcz_dispFunctionEnd(TMP_uBUOdotSvmFyUWQKMUdr);

end
