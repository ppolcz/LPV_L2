function [ret] = LPV_L2anal_LPVTools(A_fh, B_fh, C, D, K, x_lim, p_lim, dp_lim, LPVTools, varargin)
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

x_limits = x_lim';
pcz_dispFunction_num2str(x_limits)

%% LFT model

% Find out np, nw and nz
[nz,nx] = size(C);
[~,nw] = size(D);
np = size(p_lim,1);

p_cell = cell(1,np);
for i = 1:np
    name = [ 'p' num2str(i) ];
    p_cell{i} = tvreal(name, p_lim(i,:), dp_lim(i,:));
end

% Coefficient matrices of the open loop system
A_unc = A_fh(p_cell{:});
B_unc = B_fh(p_cell{:});
C_unc = plftmat(C);
D_unc = plftmat(D);

% Closed loop state transition matrix
Ak_unc = A_unc - B_unc*K;

sys_unc_LFT = ss(Ak_unc,B_unc,C_unc,D_unc);
simplify(sys_unc_LFT,'full');

%% Gridded model


p_cell = cell(1,np);
for i = 1:np
    name = [ 'p' num2str(i) ];
    p_cell{i} = pgrid(name, linspace(p_lim(i,1),p_lim(i,2),LPVTools.Resolution(i)), dp_lim(i,:));
end

rho1 = p_cell{1};
rho2 = p_cell{2};
rho3 = p_cell{3};

f1 = 1 / (rho2^2 - 2);
df1 = -2*rho2 / (rho2^2 - 2)^2;

bases = [
    basis(1,0)
%     basis(rho1,1)
%     basis(rho2,1)
%     basis(rho3,1)
%     basis(f1,df1)
    ];

% Coefficient matrices of the open loop system
A_unc = A_fh(p_cell{:});
B_unc = B_fh(p_cell{:});
C_unc = pmat(C);
D_unc = pmat(D);

% Closed loop state transition matrix
Ak_unc = A_unc - B_unc*K;

sys_unc_GRID = ss(Ak_unc,B_unc,C_unc,D_unc);

%% Computations

if LPVTools.lpvnorm_GRID
    %% -------------------------------------------------------------------------
    TMP_cwCXkgHfZmFQRzNVUlCO = pcz_dispFunctionName('LPVTools grid based', 'lpvnorm');

    [ lpv_gamma, X, info ] = lpvnorm(sys_unc_GRID, bases);
    lpvnorm_time = toc(TMP_cwCXkgHfZmFQRzNVUlCO);

    pcz_dispFunction('Solver time: <strong>%g</strong>', lpvnorm_time)
    pcz_dispFunction(2, sprintf('<strong>gamma = %g </strong>(lpvnorm)', lpv_gamma))

    pcz_dispFunctionEnd(TMP_cwCXkgHfZmFQRzNVUlCO);
    % -------------------------------------------------------------------------

    store_results('LPVTools_Results.csv', x_lim, 0, lpv_gamma, lpvnorm_time, 'tvreal, lpvnorm', 'LFT/IQC', 0)

end

if LPVTools.lpvnorm_IQC

    % -------------------------------------------------------------------------
    TMP_cwCXkgHfZmFQRzNVUlCO = pcz_dispFunctionName('LPVTools', 'lpvnorm');

    lpv_gamma = lpvnorm(sys_unc_LFT);
    lpvnorm_time = toc(TMP_cwCXkgHfZmFQRzNVUlCO);

    pcz_dispFunction('Solver time: <strong>%g</strong>', lpvnorm_time)
    pcz_dispFunction(2, sprintf('<strong>gamma = %g </strong>(lpvnorm)', lpv_gamma))

    pcz_dispFunctionEnd(TMP_cwCXkgHfZmFQRzNVUlCO);
    % -------------------------------------------------------------------------

    store_results('LPVTools_Results.csv', x_lim, 0, lpv_gamma, lpvnorm_time, 'tvreal, lpvnorm', 'LFT/IQC', 0)

end

if LPVTools.lpvwcgain_IQC

    % -------------------------------------------------------------------------
    TMP_KZWeXiYFmdpQdsgidKeG = pcz_dispFunctionName('LPVTools', 'lpvwcgain');

    wcg = lpvwcgain(sys_unc_LFT);
    lpvwcgain_time = toc(TMP_KZWeXiYFmdpQdsgidKeG);
    bounds = [wcg.LowerBound , wcg.UpperBound];

    pcz_dispFunction('Solver time: <strong>%g</strong>', lpvwcgain_time)
    pcz_dispFunction(2, 'Bounds: [<strong>%g</strong>,%g] = [%g,%g] dB ', bounds, 20*log10(bounds))

    pcz_dispFunctionEnd(TMP_KZWeXiYFmdpQdsgidKeG);
    % -------------------------------------------------------------------------

    store_results('LPVTools_Results.csv', x_lim, bounds(1), bounds(2), lpvwcgain_time, 'tvreal, lpvwcgain', 'LFT/IQC', 0)

end

%%

pcz_dispFunctionEnd(TMP_uBUOdotSvmFyUWQKMUdr);

end
