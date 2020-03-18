%% LPV_newmod_L2_gain_4D_3unc_v10_Finsler
%
%  File: LPV_newmod_L2_gain_4D_3unc_v10_Finsler.m
%  Directory: 1_PhD_projects/22_Hinf_norm
%  Author: Peter Polcz (ppolcz@gmail.com)
%
%  Created on 2018. November 11.
%

%%
% Automatically generated stuff

global SCOPE_DEPTH VERBOSE LATEX_EQNR
SCOPE_DEPTH = 0;
VERBOSE = 1;
LATEX_EQNR = 0;

try c = evalin('caller','persist'); catch; c = []; end
persist = Persist(mfilename('fullpath'), c); clear c;
persist.backup();
%clear persist

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

bounds = [
    0.1, 8,1
    0.2,12,1
    0.3,30,2
    0.4,30,3
    0.5,40,3
    0.6,50,3
    0.7,60,3
    0.8,70,3
    0.9,80,3
    1.0,90,3
    ];


x_max = bounds(6,:)';

%%

for x_max = bounds'
%%
    phi_max   = x_max(1); % radians
    v_max     = x_max(2); % meters/seconds
    omega_max = x_max(3); % radians/seconds

    x_lim = [
        -phi_max   phi_max
        -v_max     v_max
        -omega_max omega_max
        ];

    % Compute numerically the bounds for dp3
    resolution = [
        45
        45
        45
        ];
    lspace = cellfun(@(o) {linspace(o{:})+eps}, num2cell(num2cell([x_lim resolution]),2));
    grid = cell(1,nx);
    [grid{:}] = ndgrid(lspace{:});

    dx3_fh = matlabFunction(expand(simplify(F_sym(3,:)*x)));
    dx3_grid = dx3_fh(grid{:});

    p1_lim = [ sin(phi_max) / phi_max 1 ];
    p2_lim = [ cos(phi_max) 1 ];
    p3_lim = omega_max * [ -sin(phi_max) sin(phi_max) ];

    dp1_fn = @(x1,x3) ( x1*cos(x1) - sin(x1) ) / x1^2 * x3;
    dp1_lim = [ dp1_fn(phi_max,omega_max) dp1_fn(-phi_max,omega_max) ];
    dp2_lim = p3_lim;
    dp3_lim = [ floor(min(dx3_grid(:))) ceil(max(dx3_grid(:))) ];

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
    % LPV_quick_check_stability(Ak_sym, B_sym, C, D, p, p_lim)

    % Give some nice domains, where the rational matrices are well-defined and
    % on which the numerical computations will be performed (minimal generator
    % and annihilator computation).
    p_lims_comp = [
        -1 1
        -1 1
        -1 1
        ];
    pdp_lims_comp = [
        -1 1
        -1 1
        -1 1
        -1 1
        -1 1
        -1 1
        ];

    % %%

    % LPV_L2anal_Finsler(A_fh, B_fh, C, D, K, x_lim, p_expr, p_lim, dp_lim, ...
    %     p_lims_comp,pdp_lims_comp,'modelname', 'qLPV_ipend');

    % return

    % %%

%     LPV_L2anal_RCT(A_fh, B_fh, C, D, K, p_lim);

%     Use_MinLFR = 1;
%     LPV_L2anal_Masubuchi_primal(A_fh, B_fh, C, D, K, Use_MinLFR, x_lim, p_expr, p_lim, dp_lim)
%     LPV_L2anal_Masubuchi_dual(A_fh, B_fh, C, D, K, Use_MinLFR, x_lim, p_expr, p_lim, dp_lim)

    LPVTools.Resolution = [5 5 5]';
%     LPV_L2anal_Tamas_grid(A_fh, B_fh, C, D, K, x_lim, p_lim, dp_lim, LPVTools)

    %%

%     LPV_L2anal_IQCToolbox(A_fh, B_fh, C, D, K, x_lim, p_expr, p_lim, dp_lim, ...
%         p_lims_comp,pdp_lims_comp,'modelname', 'qLPV_ipend');

    LPVTools.lpvnorm_GRID = 1;
    LPVTools.lpvnorm_IQC = 0;
    LPVTools.lpvwcgain_IQC = 0;
    LPV_L2anal_LPVTools(A_fh, B_fh, C, D, K, x_lim, p_lim, dp_lim, LPVTools,'modelname', 'qLPV_ipend');

end
