%% 
%  
%  File: PAffineMatrix_self_check3_quadratic_with_affine.m
%  Directory: 7_ftools/ftools/v11/test
%  Author: Peter Polcz (ppolcz@gmail.com) 
% 
%  Created on 2018. October 23.
%  Modified on 2019. November 05. (2019a)

%%

global SCOPE_DEPTH VERBOSE LATEX_EQNR 
SCOPE_DEPTH = 0;
VERBOSE = 1;
LATEX_EQNR = 0;

% Given V(x,d). Find Pi'(x) P(d) Pi(x) = V(x,d) using P(d) as
% PAffineMatrix.

n = 2;
d_n = 3;
P_generate_symvars_v5(n,d_n);

Pi = monomials(x,1:2);
p = numel(Pi);

Theta = cellfun(@(i){ sdpvar(p,p,'symmetric') }, num2cell(0:d_n));

P = PAffineMatrix(Theta,[1;d]);

fprintf('\n\n')
pcz_info('Display channels of P(d) [they should all be sdpvariables]:\n')
P.disp_channels

P = P.to_symTheta('P');

fprintf('\n\n')
pcz_info('Display channels of P(d) after to_symTheta is called [should be symbolic symmetric matrices]:\n')
P.disp_channels

V = d1 * x1^2*x2^2 + d2 * x2^4 + d3 * x1^2 + x2^2;

% sym(P) is explicitly called, when Pi.' * P is called
V_all = expand(Pi.' * P * Pi);

c = coeffs(V_all - V, xd);

p = symvar(P.Theta);

% Nem egeszen ertem ez miert kellett annak idejen:
Theta_fh = matlabFunction(P.Theta, 'vars', {transpose(p)});

[A,B] = equationsToMatrix(c == 0, p);
A = double(A);
B = double(B);
pval = A\B;

Theta_val = Theta_fh(pval);

P.Theta = Theta_val;
display(P,'P(d) [finally computed]')

pcz_symzero(V - Pi' * P * Pi, 'V(x,d) = Pi''(x) P(d) Pi(x)')

