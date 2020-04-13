%% 
%  File: model4_randLFR.m
%  Directory: 8_published/LPV_L2/workspace
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2020. March 28. (2019b)
%

G_reset
P_init(12)

%%

RUN_ID = str2double(getenv('RUN_ID'));
if isnan(RUN_ID) || ceil(log10(RUN_ID + 1)) ~= 4
    setenv('RUN_ID', num2str(pcz_runID))
else
    setenv('RUN_ID', num2str(str2double(getenv('RUN_ID')) + 1))
end

logger = Logger(['results/' mfilename '-output.txt']);
TMP_vcUXzzrUtfOumvfgWXDd = pcz_dispFunctionName;

pcz_dispFunction2('Run ID = %s', getenv('RUN_ID'));

%%

P_generate_symvars_v10(3,1,0,0);

A0 = [
   -4.365 -0.6723 -0.3363
    7.088 -6.557  -4.601
   -2.410  7.584  -14.310
    ];

A1 = [
   -0.56081 0.85534 0.58923
    2.5333 -1.0398 -7.7373
    3.1917  1.7971 -2.5887
    ];

A2 = [
    0.66981 -1.375  -0.99093
   -2.8963  -1.5292 10.516
   -3.5777   2.8389  1.9087
    ];

B0 = [
    2.374   0.7485
    1.366   3.444
    0.9416 -9.619
    ];

B1 = [
   -0.16023 -0.35209
    0.11622 -2.4839
   -0.11058 -4.6057
    ];

B2 = [
    0.15623  0.13063
   -0.49582  4.0379
   -0.030616 0.89473
    ];

A = @(p) A0 + p*A1 + p^2*A2;

B = @(p) B0 + p*B1 + p^2*B2;

C = [
    0 1 0
    0 0 1
    ];

D = [ 0 0 ; 0 0 ];

norm_at0 = norm(ss(A(0),B(0),C,D),Inf);
norm_at1 = norm(ss(A(1),B(1),C,D),Inf);

pcz_dispFunction_scalar(norm_at0, norm_at1)

p_lim = [
    0 1
    ];

dp_lim = [
    -1 1
    ] * 1e8;

pdp_lim = [
    p_lim
    dp_lim
    ];

p_lims_comp = p_lim;
pdp_lims_comp = pdp_lim;

%%

bases = [
    1
    p1
    p1^2
    p1^3
    p1^4
    p1^5
    ];

bases_Jac = jacobian(bases,p);

modelname = sprintf('gasturbine_engine, R: %d', dp_lim(2));

%%

method1_RCT(modelname,A,B,C,D,p_lim)

method0_grid_LPVTools(modelname,A,B,C,D,p_lim,dp_lim,bases,bases_Jac,[5 5 5]);

% Greedy grid 
method0_grid_Wu1995(modelname,A,B,C,D,p_lim,dp_lim,bases,bases_Jac,'res_max',5');
method0_grid_Wu1995(modelname,A,B,C,D,p_lim,dp_lim,bases,bases_Jac,'res_max',15');

% As proposed by Wu (1995,1996)
% method0_grid_Wu1995(modelname,A,B,C,D,p_lim,dp_lim,bases,bases_Jac,'res_max',1000,'delta',1e-4,'T',10000);
% method0_grid_Wu1995(modelname,A,B,C,D,p_lim,dp_lim,bases,bases_Jac,'res_max',1000,'delta',1e-3,'T',1000);
% method0_grid_Wu1995(modelname,A,B,C,D,p_lim,dp_lim,bases,bases_Jac,'res_max',1000,'delta',1e-2,'T',100);
% method0_grid_Wu1995(modelname,A,B,C,D,p_lim,dp_lim,bases,bases_Jac,'res_max',1000,'delta',1e-1,'T',100);
% method0_grid_Wu1995(modelname,A,B,C,D,p_lim,dp_lim,bases,bases_Jac,'res_max',1000,'delta',1e-1,'T',10);

method2_descriptor_primal(modelname,A,B,C,D,p_lim,dp_lim)
method2_descriptor_dual(modelname,A,B,C,D,p_lim,dp_lim,1)
method2_descriptor_dual(modelname,A,B,C,D,p_lim,dp_lim,0)

% IQC/LFT approaches for LPV with rate-bounded parameters
method3_IQC_LFT_IQCToolbox(modelname,A,B,C,D,p_lim,dp_lim);
method3_IQC_LFT_LPVTools(modelname,A,B,C,D,p_lim,dp_lim);

method4_authors_old_symbolical(modelname,A,B,C,D,p_lim,dp_lim);

% Imported variables to the base workspace: Q, dQ, PI_x, gamma
method5_proposed_approach(modelname,A,B,C,D,p_lim,dp_lim,p_lim,pdp_lim,[]);

%%

pcz_dispFunctionEnd(TMP_vcUXzzrUtfOumvfgWXDd);
logger.stoplog
