%% PGenAffineMatrix_basic_tests
%  
%  File: PGenAffineMatrix_basic_tests.m
%  Directory: 7_ftools/ftools/v11/test
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2019. September 06. (2019a)
%

%%

P_generate_symvars_v5(3,0);

Theta_ = round(10*rand(4,6));

P = PGenAffineMatrix(Theta_,[x2 sin(x1) ],x)

P(1,2,3)

%%

Theta = sdpvar([3 3]);

P = PGenAffineMatrix(Theta,[x2,sin(x1)],x)

P(0,2,3)