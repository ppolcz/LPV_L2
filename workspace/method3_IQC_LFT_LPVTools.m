function [ret] = method3_IQC_LFT_LPVTools(A_fh, B_fh, C, D, x_lim, p_lim, dp_lim, opts, varargin)
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

sys_unc_LFT = ss(A_unc,B_unc,C_unc,D_unc);
simplify(sys_unc_LFT,'full');

%% Gridded model


p_cell = cell(1,np);
for i = 1:np
    name = [ 'p' num2str(i) ];
    p_cell{i} = pgrid(name, linspace(p_lim(i,1),p_lim(i,2),opts.Resolution(i)), dp_lim(i,:));
end

r1 = p_cell{1};
r2 = p_cell{2};
r3 = p_cell{3};

% bases = [
%                    1
%                   p1
%                   p2
%                   p3
%  -(p1*p2)/(p2^2 - 7)
%       -p1/(p2^2 - 7)
%       -p2/(p2^2 - 7)
%      p2^2/(p2^2 - 7)
%  -(p2*p3)/(p2^2 - 7)
%       -p3/(p2^2 - 7)
%     ];

bases = [
    basis(                   1, 'p1',              0, 'p2',                                         0, 'p3',              0)
    basis(                  r1, 'p1',              1, 'p2',                                         0, 'p3',              0)
    basis(                  r2, 'p1',              0, 'p2',                                         1, 'p3',              0)
    basis(                  r3, 'p1',              0, 'p2',                                         0, 'p3',              1)
    basis( -(r1*r2)/(r2^2 - 7), 'p1', -r2/(r2^2 - 7), 'p2',  (2*r1*r2^2)/(r2^2 - 7)^2 - r1/(r2^2 - 7), 'p3',              0)
    basis(      -r1/(r2^2 - 7), 'p1',  -1/(r2^2 - 7), 'p2',                    (2*r1*r2)/(r2^2 - 7)^2, 'p3',              0)
    basis(      -r2/(r2^2 - 7), 'p1',              0, 'p2',      (2*r2^2)/(r2^2 - 7)^2 - 1/(r2^2 - 7), 'p3',              0)
    basis(     r2^2/(r2^2 - 7), 'p1',              0, 'p2', (2*r2)/(r2^2 - 7) - (2*r2^3)/(r2^2 - 7)^2, 'p3',              0)
    basis( -(r2*r3)/(r2^2 - 7), 'p1',              0, 'p2',  (2*r2^2*r3)/(r2^2 - 7)^2 - r3/(r2^2 - 7), 'p3', -r2/(r2^2 - 7))
    basis(      -r3/(r2^2 - 7), 'p1',              0, 'p2',                    (2*r2*r3)/(r2^2 - 7)^2, 'p3',  -1/(r2^2 - 7))
    ];

% Coefficient matrices of the open loop system
A_unc = A_fh(p_cell{:});
B_unc = B_fh(p_cell{:});
C_unc = pmat(C);
D_unc = pmat(D);

sys_unc_GRID = ss(A_unc,B_unc,C_unc,D_unc);

%% Computations

if opts.lpvnorm_GRID
    %% -------------------------------------------------------------------------
    TMP_cwCXkgHfZmFQRzNVUlCO = pcz_dispFunctionName('LPVTools grid based', 'lpvnorm');

    [ lpv_gamma, X, info ] = lpvnorm(sys_unc_GRID, bases);
    lpvnorm_time = toc(TMP_cwCXkgHfZmFQRzNVUlCO);

    pcz_dispFunction('Solver time: <strong>%g</strong>', lpvnorm_time)
    pcz_dispFunction(2, sprintf('<strong>gamma = %g </strong>(lpvnorm)', lpv_gamma))

    pcz_dispFunctionEnd(TMP_cwCXkgHfZmFQRzNVUlCO);
    % -------------------------------------------------------------------------

    store_results('LPVTools_Results.csv', x_lim, 0, lpv_gamma, lpvnorm_time, 'tvreal, lpvnorm', 'grid', 0)
    store_results('Results_All.csv', x_lim, 0, lpv_gamma, lpvnorm_time, 'tvreal, lpvnorm', 'LPVTools - grid', 0)

end

if opts.lpvnorm_IQC

    % -------------------------------------------------------------------------
    TMP_cwCXkgHfZmFQRzNVUlCO = pcz_dispFunctionName('LPVTools LFT/IQC', 'lpvnorm');

    lpv_gamma = lpvnorm(sys_unc_LFT);
    lpvnorm_time = toc(TMP_cwCXkgHfZmFQRzNVUlCO);

    pcz_dispFunction('Solver time: <strong>%g</strong>', lpvnorm_time)
    pcz_dispFunction(2, sprintf('<strong>gamma = %g </strong>(lpvnorm)', lpv_gamma))

    pcz_dispFunctionEnd(TMP_cwCXkgHfZmFQRzNVUlCO);
    % -------------------------------------------------------------------------

    store_results('LPVTools_Results.csv', x_lim, 0, lpv_gamma, lpvnorm_time, 'tvreal, lpvnorm', 'LFT/IQC', 0)
    store_results('Results_All.csv', x_lim, 0, lpv_gamma, lpvnorm_time, 'tvreal, lpvnorm', 'LPVTools - LFT/IQC', 0)

end

if opts.lpvwcgain_IQC

    % -------------------------------------------------------------------------
    TMP_KZWeXiYFmdpQdsgidKeG = pcz_dispFunctionName('LPVTools LFT/IQC', 'lpvwcgain');

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

    try    
        store_results('LPVTools_Results.csv', x_lim, bounds(1), bounds(2), lpvwcgain_time, 'tvreal, lpvwcgain', 'LFT/IQC', 0)
        store_results('Results_All.csv', x_lim, bounds(1), bounds(2), lpvwcgain_time, 'tvreal, lpvwcgain', 'LPVTools - LFT/IQC', 0)
    catch ex
        getReport(ex)
    end

end

%%

pcz_dispFunctionEnd(TMP_uBUOdotSvmFyUWQKMUdr);

end
