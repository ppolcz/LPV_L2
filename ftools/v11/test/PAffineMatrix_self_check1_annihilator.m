%% 
%  File: PAffineMatrix_self_check1_annihilator.m
%  Directory: 1_PhD_projects/00_my_toolboxes/FinslerTools/v11/test
%  Author: Peter Polcz (ppolcz@gmail.com) 
% 
%  Created on 2018. October 23.

%% 

global SCOPE_DEPTH VERBOSE LATEX_EQNR 
SCOPE_DEPTH = 0;
VERBOSE = 1;
LATEX_EQNR = 0;

TMP_VDFhUxkUNToJIgvoLoUQ = pcz_dispFunctionName;

n = 5;
P_generate_symvars_v5(n,0);

Theta = rand(4,2*n);
Theta = ( Theta < 0.07 ) - (Theta > 0.93);

N = PAffineMatrix(Theta,x,[x2;x4])
N = N.generate_symbolic
N(rand(2,1))

N.disp_channels

pcz_dispFunctionEnd(TMP_VDFhUxkUNToJIgvoLoUQ);
