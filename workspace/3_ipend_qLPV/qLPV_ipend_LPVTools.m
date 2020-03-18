%% qLPV_ipend_LPVTools
%  
%  File: qLPV_ipend_LPVTools.m
%  Directory: 1_PhD_projects/22_Hinf_norm/qLPV_ipend
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2020. March 10. (2019b)
%

%%

global SCOPE_DEPTH VERBOSE LATEX_EQNR
SCOPE_DEPTH = 0;
VERBOSE = 1;
LATEX_EQNR = 0;

%%

O = @(varargin) zeros(varargin{:});
I = @(varargin) eye(varargin{:});
He = he;

%% Nonlinear model (as a starting point)

P_generate_symvars_v10(3,3,1,1);

% Inverted pendulum model and parameters taken from [Szederkenyi, Hangos,
% Bokor and Vamos, 2002, Linear output selection for feedback
% linearization]
m = 0.5; % [kg], mass of the rod
M = 0.5; % [kg], mass of the car
l = 1;   % [m], length of the rod
g = 10;  % [m/s^2], gravitational acceleration

% Constants describing the motor characteristics:
alpha = 0.4;
beta = 0.3;

% Nonlinear model:
% dx = f(x) + G(x)*u, where f(x) = F(x)*x
F_sym = [
    0                                                0                                      1
    -m*g*sin(x1)*cos(x1) / ( x1*(M+m*sin(x1)^2) )   -beta / ( M+m*sin(x1)^2 )               m*l*x3*sin(x1) / ( M+m*sin(x1)^2 )
    (M+m)*g*sin(x1) / ( x1*l*(M+m*sin(x1)^2) )       beta*cos(x1) / ( l*(M+m*sin(x1)^2) )  -m*x3*sin(x1)*cos(x1) / ( M+m*sin(x1)^2 )
    ];
G_sym = [
    0
    alpha / ( M+m*sin(x1)^2 )
    alpha*cos(x1) / ( l*(M+m*sin(x1)^2) )
    ];

% Linearization around x = 0:
A = double(subs(jacobian(F_sym*x,x),x,0*x));
B = double(subs(G_sym,x,0*x));

% Design a static state feedback:
K = lqr(A,B,eye(3),1);

% State equation of the closed loop system:
f_cls = simplify(F_sym*x) - G_sym*K*x;

% Quick check using a simple simulation:
%{
%%

f_ode = matlabFunction(f_cls, 'vars', {t x});
G_ode = matlabFunction(G_sym,'vars',{t,x});
% ode45(f_ode, [0 10], [1.1 0 0]');

u = @(t) interp1([0,1-eps,1+eps,1000],[1,1,0,0],t,'linear');
f_ode_disturbance = @(t,x) f_ode(t,x) + G_ode(t,x)*u(t);
[tt,xx] = ode45(f_ode_disturbance, [0,10],[0 0 0]');
yy = xx(:,1);

L2_NORM_u = 1;
L2_NORM_y = sum(yy(1:end-1).^2 .* diff(tt));

plot(tt,xx)

%%
%}

%% qLPV model of the inverted pendulum balance system

% Build up a qLPV model using artificial state-dependent parameters:
p_expr = [
    sin(x1)/x1  % = p1
    cos(x1)     % = p2
    x3*sin(x1)  % = p3
    ];

%%% 
% Coefficient matrices of the qLPV model.
% 
%  dx = A(p)*x + B(p)*u
%   y = C*x
%   u = -K*x + w (disturbance)
%  
A_fh = @(p1,p2,p3) [
    0                                 0                              1
    -m*g*p1*p2 / ( M+m-m*p2^2 )      -beta / ( M+m-m*p2^2 )          m*l*p3 / ( M+m-m*p2^2 )
    (M+m)*g*p1 / ( l*(M+m-m*p2^2) )   beta*p2 / ( l*(M+m-m*p2^2) )  -m*p2*p3 / ( M+m-m*p2^2 )
    ];

B_fh = @(p1,p2,p3) [
    0
    alpha / ( M+m-m*p2^2 )
    alpha*p2 / ( M+m-m*p2^2 )
    ];

C = [ 1 0 0 ];

D = 0;

% Compute the symbolical expressions of the coefficient matrices to check
% correctness.
A_sym = A_fh(p1,p2,p3);
B_sym = B_fh(p1,p2,p3);
Ak_sym = A_sym - B_sym*K;

% Check symbolically that the qLPV model with the given parameter
% expressions corresponds to the nonlinear model.
pcz_symeq_report(F_sym, subs(A_sym,p,p_expr), 'F(x) == A(p=p_expr)')
pcz_symeq_report(G_sym, subs(B_sym,p,p_expr), 'G(x) == B(p=p_expr)')

%% Define parameter bounds

phi_max = 0.2; % radians
v_max = 10; % meters/seconds
omega_max = 30; % radians/seconds

p1_lim = [ sin(phi_max) / phi_max 1 ];
p2_lim = [ cos(phi_max) 1 ];
p3_lim = omega_max * [ -sin(phi_max) sin(phi_max) ];

dp1_fn = @(x1,x3) ( x1*cos(x1) - sin(x1) ) / x1^2 * x3;
dp1_lim = [ dp1_fn(phi_max,omega_max) dp1_fn(-phi_max,omega_max) ];
dp2_lim = p3_lim;
dp3_lim = [-234 234];

x_lim = [
    -phi_max phi_max
    -v_max v_max
    -omega_max omega_max
    ];
    

p_lim = [
    p1_lim
    p2_lim
    p3_lim
    ];

dp_lim = [
    dp1_lim
    dp2_lim
    dp3_lim
    ];

pdp_lim = [
    p_lim
    dp_lim
    ];

% Check stability in random point
LPV_quick_check_stability(Ak_sym, B_sym, C, D, p, p_lim)

%%

F_fh = @(p1,p2,p3) [
    A_fh(p1,p2,p3) , B_fh(p1,p2,p3)
    C              , D
    ];

F_lfr = F_fh(p_cell{:});
% F_lfr = minlfr(F_lfr);

%%

pnum_cell = num2cell(rand(1,np));
F_rand = F_fh(pnum_cell{:});
[nxz,nxw] = size(F_rand);

nz = nxz - nx;
nw = nxw - nx;

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

sys_unc = ss(Ak_unc,B_unc,C_unc,D_unc);
simplify(sys_unc,'full');


% -------------------------------------------------------------------------
TMP_cwCXkgHfZmFQRzNVUlCO = pcz_dispFunctionName('LPVTools', 'lpvnorm');

lpv_gamma = lpvnorm(sys_unc);
pcz_dispFunction('Solver time: <strong>%g</strong>', toc(TMP_cwCXkgHfZmFQRzNVUlCO))
pcz_dispFunction(2, sprintf('<strong>gamma = %g </strong>(lpvnorm)', lpv_gamma))

pcz_dispFunctionEnd(TMP_cwCXkgHfZmFQRzNVUlCO);
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
TMP_KZWeXiYFmdpQdsgidKeG = pcz_dispFunctionName('LPVTools', 'lpvwcgain');

wcg = lpvwcgain(sys_unc);
bounds = [wcg.LowerBound , wcg.UpperBound];

pcz_dispFunction('Solver time: <strong>%g</strong>', toc(TMP_KZWeXiYFmdpQdsgidKeG))
pcz_dispFunction(2, 'Bounds: [<strong>%g</strong>,%g] = [%g,%g] dB ', bounds, 20*log10(bounds))

pcz_dispFunctionEnd(TMP_KZWeXiYFmdpQdsgidKeG);
% -------------------------------------------------------------------------
