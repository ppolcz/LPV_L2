%% 
%
%  File: PAffineMatrix_self_check2_quadratic.m
%  Directory: 1_PhD_projects/00_my_toolboxes/FinslerTools/v11/test
%  Author: Peter Polcz (ppolcz@gmail.com) 
% 
%  Created on 2018. October 23.

%%

global SCOPE_DEPTH VERBOSE LATEX_EQNR 
SCOPE_DEPTH = 0;
VERBOSE = 1;
LATEX_EQNR = 0;

TMP_SYDZCArwzdUdwQuFQUWj = pcz_dispFunctionName;

n = 5;
P_generate_symvars_v5(n,0);

Theta = sdpvar(4,2*n);

P = PAffineMatrix(Theta,x,[x2;x4],'name','P');
P = P.to_symTheta('P');

pcz_dispFunction2('In the followings you should see the symbolical values of the channels.')
pcz_dispFunction2(evalc('P.disp_channels'))
pcz_dispFunction2(evalc('display(P,''P(x)'')'))

ujP = P * [1 2 ; 2 0];
pcz_symzero(P.Sym * [1 2 ; 2 0] - ujP.Sym, 'P(x)*[1 2;2 0] = ujP(x)')

P_sym = sym(P);

pcz_symzero(subs(P_sym,[x1;x2;x3;x4;x5],[0;1;0;2;0]) - P([1;2]), 'P(val) = subs(sym(P),x,val)')

pcz_dispFunctionEnd(TMP_SYDZCArwzdUdwQuFQUWj);

%%
