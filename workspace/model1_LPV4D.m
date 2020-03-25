%% LPV4D_main
%  
%  File: LPV4D_main.m
%  Directory: workspace/1_comp_LPV
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2020. March 20. (2019b)
%

global SCOPE_DEPTH VERBOSE
SCOPE_DEPTH = 0;
VERBOSE = 1;

logger = Logger('results/qLPV_ipend_main-output.txt');
TMP_vcUXzzrUtfOumvfgWXDd = pcz_dispFunctionName;

RUN_ID = str2double(getenv('RUN_ID'));
if isnan(RUN_ID) || ceil(log10(RUN_ID + 1)) ~= 4
    setenv('RUN_ID', num2str(pcz_runID))
else
    setenv('RUN_ID', num2str(str2double(getenv('RUN_ID')) + 1))
end

%%

p_lim = [
    -1 2
    -1 2
    0 2
    ];

dp_lim = [
    -10 10
    -1 1
    -5 5
    ];

[p_lfr,p_lfr_cell] = pcz_generateLFRStateVector('p',p_lim);

A_fh = @(p1,p2,p3) [
    -3+p1       3+p1/(p3^2 + 0.5*p1 + 1)     0.1*p3   p3^2*p2
    0          -1-p2^2                       5        0
    -1/(5-p2)  0                             -4+p1    0
    0          0.1                           0        -5+1/(p1 + 2)
    ];
    
C_fh = @(p1,p2,p3) [
    1/(5-p2)   0                             0        0
    0          0                             p1+1     0
    ];

B_fh = @(p1,p2,p3) [
    0      0
    1+p2^2 0
    0      0
    0      2+p1/(p1 + 2)
    ];

D = [
    0      0
    0      0
    ];

M_lfr = minlfr([
    A_fh(p_lfr_cell{:}) B_fh(p_lfr_cell{:})
    C_fh(p_lfr_cell{:}) D
    ])

%%

% LPV_quick_check_stability(A_fh, B_fh, C_fh, D, p_lim)

p_lims_comp = p_lim;

pdp_lims_comp = [
    p_lims_comp
    dp_lim
    ];
    
%%

% LPV_L2anal_Tamas_grid(A_fh, B_fh, C_fh, D, p_lim, dp_lim, [5 5 5])
% 
% LPV_L2anal_RCT(A_fh, B_fh, C_fh, D, p_lim);
% 
% LPV_L2anal_IQCToolbox(A_fh, B_fh, C_fh, D, p_lim, dp_lim);

% LPV_L2anal_Finsler(A_fh, B_fh, C_fh, D, p_lim, dp_lim, p_lims_comp,pdp_lims_comp);
% LPV_L2anal_Finsler_old_symbolical(A_fh, B_fh, C_fh, D, p_lim, dp_lim, p_lims_comp, pdp_lims_comp)

% Use_MinLFR = 1;
% LPV_L2anal_Masubuchi_primal(A_fh, B_fh, C_fh, D, Use_MinLFR, p_lim, dp_lim)
% LPV_L2anal_Masubuchi_dual(A_fh, B_fh, C_fh, D, Use_MinLFR, p_lim, dp_lim)

% LPVTools.lpvnorm_GRID = 1;
% LPVTools.Resolution = [5 5 5]';
% LPVTools.lpvnorm_IQC = 1;
% LPVTools.lpvwcgain_IQC = 1;
% LPV_L2anal_LPVTools(A_fh, B_fh, C_fh, D, p_lim, dp_lim, LPVTools);

pcz_dispFunctionEnd(TMP_vcUXzzrUtfOumvfgWXDd);
logger.stoplog
