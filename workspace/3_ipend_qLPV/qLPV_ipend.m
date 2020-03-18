%% 
%  File: qLPC_ipend.m
%  Directory: 1_PhD_projects/22_Hinf_norm/qLPV_ipend
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2020. March 07. (2019b)
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

phi_max = 0.3; % radians
v_max = 10; % meters/seconds
omega_max = 30; % radians/seconds

p1_lim = [ sin(phi_max) / phi_max 1 ];
p2_lim = [ cos(phi_max) 1 ];
p3_lim = omega_max * [ -sin(phi_max) sin(phi_max) ];

dp1_fn = @(x1,x3) ( x1*cos(x1) - sin(x1) ) / x1^2 * x3;
dp1_lim = [ dp1_fn(phi_max,omega_max) dp1_fn(-phi_max,omega_max) ];
dp2_lim = p3_lim;

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
    ];

pdp_lim = [
    p_lim
    dp_lim
    ];

% Check stability in random point
LPV_quick_check_stability(Ak_sym, B_sym, C, D, p, p_lim)

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
    ];

%%
% We used LFR Toolbox objects (class `lfr'), which are often wrapped into a
% class `plfr' object.

% Generate symbolic lfr variables
[x_lfr,x_cell] = pcz_generateLFRStateVector('x',nx);
[u_lfr,u_cell] = pcz_generateLFRStateVector('u',nw);
[p_lfr,p_cell] = pcz_generateLFRStateVector('p',p_lim);
[dp_lfr,dp_cell] = pcz_generateLFRStateVector('dp',dp_lim);

dp_lfr_all = [ dp_lfr ; 0 ];
dp_cell_all = [ dp_cell { 0 } ];

p_sym = sym(plfr(p_lfr));
dp_sym = sym(plfr(dp_lfr));
dp_sym_all = sym(plfr(dp_lfr_all));
pdp_sym = [ p_sym ; dp_sym ];


% LFR realization of matrix B(p):
B_lfr = B_fh(p_cell{:});

% LFR realization of matrix A(p) - B(p)*K:
Ak_lfr = A_fh(p_cell{:}) - B_lfr*K;

% -----------------------

% Model matrices Fij collected in cell F{:,:}
F = cell(4,4);

%%%
% LFR of matrix A(p) = (F11 F12)*PI_x, where PI_x = [ eye(3) ; PI_1 ]
% 
% Subscript |_x| is not used in the manuscript. |PI_x|, |S_x| and |iS_x| in
% the script corresponds to $\Pi$, $S$ and $S^{-1}$ in the manuscript.
Ak_plfr = plfr(Ak_lfr);
[F{1,1},F{1,3},F{3,1},F{3,3},PI_x,PI_1] = deal(Ak_plfr.A, Ak_plfr.B, Ak_plfr.C, Ak_plfr.D, Ak_plfr.generatePI, Ak_plfr.generatePI1);
pcz_symzero_report(Ak_sym - [ F{1,1} F{1,3} ] * sym(PI_x),sprintf('Generator form (m = %d) of A(p)-B(p)*K is OK', PI_1.ny));

% Minimal generator for PI_x
[S_x,PI_x,iS_x,Ker_x] = P_mingen_for_LFR(PI_x,'lims',p_lims_comp);
[F{1,1},F{1,3},F{3,1},F{3,3}] = pcz_split_matrix([F{1,1} F{1,3} ; F{3,1} F{3,3}] * S_x, [nx Inf], [nx Inf]);
m1 = size(F{1,3},2);
mx = nx + m1;
pcz_symzero_report(Ak_sym - [ F{1,1} F{1,3} ] * sym(PI_x),'Minimal generator form of A(p)-B(p)*K is OK');

% Reload PI_1 corresponding to the minimal generator PI_x.
PI_1 = plfr(iS_x(nx+1:end,nx+1:end) * PI_1); % plfr object

% C(x,p)
F{2,1} = C;
F{2,3} = O(nz,m1);

% -----------------------

% B(x,p)
B_plfr = plfr(B_lfr);
[F{1,2},F{1,4},F{4,2},F{4,4},PI_u,PI_2] = deal(B_plfr.A, B_plfr.B, B_plfr.C, B_plfr.D, B_plfr.generatePI,B_plfr.generatePI1);
pcz_symzero_report(B_sym - [ F{1,2} F{1,4} ] * sym(PI_u),'Generator form of B(p) is OK');

% Minimal generator for PI_u
[S_u,PI_u,iS_u,Ker_u] = P_mingen_for_LFR(PI_u,'lims',[-1 1]);
[F{1,2},F{1,4},F{4,2},F{4,4}] = pcz_split_matrix([F{1,2},F{1,4};F{4,2},F{4,4}] * S_u, [nx Inf], [nw Inf]);
m2 = size(F{1,4},2);
mu = nw + m2;
pcz_symzero_report(B_sym - [ F{1,2} F{1,4} ] * sym(PI_u),'Minimal generator form of B(p) is OK');

% Reload PI_2 corresponding to the minimal generator PI_u.
PI_2 = plfr(iS_u(nw+1:end,nw+1:end) * PI_2); % plfr object

% D(x,p)
F{2,2} = D;
F{2,4} = O(nz,m2);

%% Matrices of the dissipativity inequality

% Upper-left block matrix of $\Pi_a$ (LFR Toolbox object)
PI_ax = [ % lfr object
    PI_x.lfrtbx_obj
    PI_1 * ( F{1,1} + F{1,3}*PI_1 )
    PI_1.diff(p_cell,dp_cell_all)
    ];

% Lower-right block matrix of $\Pi_a$ (LFR Toolbox object)
PI_au = [ % lfr object
    PI_u.lfrtbx_obj
    PI_1 * ( F{1,2} + F{1,4}*PI_2 )
    ];

PI_a = blkdiag(PI_ax,PI_au);

% Minimal generator for PI_a
PI_a_old = plfr(PI_a,[p_lfr;dp_lfr]);
[S_a,PI_a,iS_a,Ker_a] = P_mingen_for_LFR(PI_a_old,'lims',pdp_lims_comp);
pcz_fhzero_report(@(pdp) PI_a_old(pdp)-S_a*PI_a(pdp), pdp_sym);

Aa = [
    F{1,1}   F{1,3}   O(nx,m1) O(nx,m1)
    O(m1,nx) O(m1,m1) I(m1)    I(m1)   
    ];

Ca = [ 
    F{2,1}   F{2,3}   O(nz,m1) O(nz,m1)
    ];

Ba = [
    F{1,2}   F{1,4}   O(nx,m1)
    O(m1,nw) O(m1,m2) I(m1)
    ];

Da = [ 
    F{2,2}   F{2,4}   O(nz,m1) 
    ];

Ea = [ I(nx+m1) O(nx+m1,m1+m1) ];

Ga = [ I(nw) O(nw,m2+m1) ];
    
%% Construct annihilators

% Find `strategically good' points:
%{
best.minerr = 5.0759e-15;
for i = 1:100
    [N,N_samples,N_sampleserr,N_err] = P_affine_annihilator_for_LFR(PI_x,p_lfr,'sym',1,'lims',p_lims_comp);
    
    if N_err < best.minerr
        best.minerr = N_err;
        best.N_samples = N_samples;
        best.N_sampleserr = N_sampleserr;
        best.N = N;
    end
end

best.minerra = 2.9667e-12;
for i = 1:100
    [Na,Na_samples,Na_sampleserr,Na_err] = P_affine_annihilator_for_LFR(PI_a,[p_lfr;dp_lfr],'lims',pdp_lims_comp);
    
    if Na_err < best.minerra
        best.minerra = Na_err
        best.Na_samples = Na_samples;
        best.Na_sampleserr = Na_sampleserr;
        best.Na = Na;
    end
end
%}

% Generate annihilator N(p) for PI(p).
% N_samples = [
%       -0.308   -0.782   -0.661    0.867   -0.956    0.628   -0.886   -0.811   -0.284    0.284    0.071   -0.845   -0.477   -0.510   -0.481   -0.474   -0.649    0.819
%       -0.775   -0.533    0.721    0.649   -0.091    0.533   -0.285   -0.520   -0.577   -0.107   -0.639   -0.695   -0.734    0.698    0.908   -0.497   -0.970    0.263
%        0.238   -0.106    0.452    0.811    0.737   -0.718   -0.100    0.157   -0.726   -0.071    0.406    0.469   -0.988   -0.798    0.955    0.987    0.462   -0.508
%     ];
% N_sampleserr = [
%        0.217   -0.192    0.186    0.978   -0.215    0.237   -0.938   -0.123    0.470   -0.932   -0.407   -0.664    0.597   -0.815   -0.534    0.352   -0.660    0.850
%       -0.259   -0.529   -0.415   -0.396    0.602    0.339   -0.427   -0.561    0.068   -0.607   -0.821   -0.465   -0.165   -0.484    0.549   -0.609    0.244    0.364
%        0.600   -0.351    0.587   -0.004   -0.083   -0.692   -0.030   -0.519   -0.674    0.384    0.193   -0.539    0.404   -0.112    0.699   -0.265   -0.970   -0.498
%     ];
% [N,~,~,N_err] = P_affine_annihilator_for_LFR(PI_x,p_lfr,'sym',1,...
%     'samples',N_samples,'sampleserr',N_sampleserr);
[N,~,~,N_err] = P_affine_annihilator_for_LFR(PI_x,p_lfr,'sym',1,'lims',p_lims_comp);
pcz_fhzero_report(@(p) N(p)*PI_x(p), p_sym, N_err, 'Annihilator Na');


% Generate annihilator Na(p,dp) for PI_a(p,dp)
% Na_samples = [
%        0.810    0.769    0.725    0.831    0.461   -0.307    0.139    0.243   -0.921    0.917    0.078   -0.819   -0.696    0.167   -0.005    0.210   -0.904    0.050   -0.752   -0.489   -0.715    0.455    0.751   -0.112    0.804   -0.863    0.131    0.959   -0.820   -0.875    0.994    0.841    0.764   -0.664    0.840    0.443   -0.440    0.124   -0.106    0.862   -0.123   -0.027    0.370   -0.845   -0.827    0.561   -0.553   -0.002    0.381    0.408   -0.451   -0.950   -0.317   -0.584   -0.496   -0.923    0.492
%       -0.240    0.165    0.744   -0.368    0.240    0.875    0.078    0.451    0.018    0.786   -0.058   -0.609    0.624   -0.768   -0.930    0.525    0.281    0.490   -0.154    0.923   -0.170    0.310   -0.392    0.784    0.661    0.062   -0.975   -0.762   -0.527    0.238    0.529   -0.275   -0.569   -0.567    0.929   -0.584    0.880   -0.921    0.863    0.655   -0.854    0.079   -0.917    0.489    0.514   -0.447    0.152   -0.265   -0.920   -0.268    0.999    0.556    0.492   -0.870   -0.355    0.506   -0.474
%        0.580    0.662    0.578    0.379   -0.886    0.285    0.813    0.996   -0.681    0.569   -0.645   -0.341   -0.859    0.505    0.103   -0.410    0.322   -0.923    0.941   -0.376    0.317    0.592   -0.387   -0.089   -0.015    0.370    0.209   -0.952   -0.757    0.573    0.013    0.709    0.332   -0.627   -0.758    0.413   -0.113   -0.210   -0.352    0.464   -0.328    0.552    0.066    0.486    0.030   -0.522   -0.023   -0.832    0.240    0.620    0.991    0.190    0.088    0.019   -0.084   -0.381    0.369
%       -0.881   -0.417   -0.885   -0.867    0.284    0.030    0.305   -0.798   -0.850    0.448   -0.206   -0.424    0.772   -0.444    0.441   -0.094    0.374    0.200   -0.415   -0.140    0.462    0.850   -0.877   -0.561    0.583   -0.807   -0.875    0.792    0.988    0.059   -0.904    0.071   -0.459   -0.573    0.659   -0.900    0.827   -0.338    0.949   -0.959    0.357    0.269   -0.509   -0.647    0.814   -0.661    0.373    0.862    0.037   -0.849   -0.541   -0.810    0.737    0.874   -0.514    0.420   -0.213
%       -0.204   -0.177    0.278   -0.380   -0.156   -0.226   -0.631   -0.597   -0.079    0.766    0.845   -0.976   -0.815    0.596    0.288    0.292    0.483   -0.109   -0.643    0.596    0.692   -0.949   -0.379    0.573    0.412    0.528   -0.121   -0.216    0.142   -0.854    0.074    0.118    0.586   -0.252    0.828   -0.140   -0.262   -0.538    0.460   -0.667   -0.532    0.721    0.714   -0.187   -0.378    0.872    0.472   -0.615   -0.640    0.660   -0.844    0.405    0.440   -0.721   -0.487   -0.223    0.368
%     ];
% Na_sampleserr = [
%        0.969    0.132    0.922   -0.349    0.113    0.925   -0.332    0.295    0.209   -0.046    0.293   -0.995    0.876   -0.380    0.337    0.677    0.087    0.966    0.098   -0.956    0.259   -0.254    0.519   -0.163   -0.427   -0.737    0.436    0.033   -0.529   -0.296   -0.962   -0.988    0.431   -0.943   -0.953    0.087    0.794    0.611   -0.639   -0.250   -0.222    0.373   -0.896   -0.682    0.417    0.800   -0.997   -0.559    0.848   -0.868   -0.646    0.186   -0.245    0.308   -0.210   -0.128    0.761
%        0.569    0.605    0.006   -0.097   -0.821    0.004   -0.397   -0.856    0.643    0.178   -0.946    0.602   -0.131    0.214    0.254    0.227    0.565    0.276   -0.831   -0.939   -0.849    0.834    0.251    0.754    0.456   -0.900    0.375    0.164    0.791   -0.620   -0.637    0.913    0.255    0.034   -0.052    0.417   -0.203   -0.033   -0.148   -0.125    0.883    0.444   -0.147    0.324   -0.100   -0.671    0.106   -0.228    0.793   -0.113    0.921   -0.071    0.059    0.059    0.495   -0.025   -0.740
%        0.310    0.520   -0.590   -0.141    0.120    0.287   -0.210   -0.259   -0.827   -0.009    0.794   -0.954   -0.944   -0.446   -0.387    0.163    0.617    0.461    0.165    0.583   -0.013   -0.628   -0.058   -0.550    0.093    0.056    0.234    0.255    0.687   -0.271    0.868    0.068   -0.707   -0.649   -0.717    0.165    0.474   -0.186    0.376    0.754   -0.848    0.947    0.917    0.795    0.932    0.912   -0.366   -0.639   -0.813   -0.775    0.767   -0.304    0.987   -0.967    0.360    0.321   -0.823
%       -0.037    0.106    0.860   -0.532    0.695   -0.694   -0.675    0.441    0.798    0.753   -0.397   -0.499   -0.278   -0.393   -0.699   -0.167   -0.331   -0.696    0.790    0.499   -0.917    0.754   -0.408   -0.367    0.982    0.781    0.115   -0.217   -0.505    0.830   -0.746   -0.491   -0.719   -0.107    0.895   -0.906    0.297    0.727    0.772    0.911   -0.816    0.874    0.784    0.620   -0.191   -0.348    0.490   -0.192   -0.316    0.556    0.396    0.015   -0.393    0.024    0.003    0.061   -0.047
%       -0.675   -0.835   -0.892   -0.611    0.317    0.727    0.766   -0.189   -0.625    0.669    0.160   -0.990   -0.594   -0.247    0.242    0.365    0.768   -0.233   -0.188   -0.844    0.729    0.349    0.596   -0.437   -0.948   -0.950   -0.932    0.174    0.647   -0.597   -0.473   -0.468   -0.465    0.365   -0.506    0.104    0.333   -0.489    0.114    0.845    0.004    0.079   -0.414   -0.975    0.451   -0.060    0.343    0.114    0.098    0.508    0.174    0.035    0.317   -0.298    0.431    0.604   -0.057
%     ];
% [Na,~,~,Na_err] = P_affine_annihilator_for_LFR(PI_a,[p_lfr;dp_lfr],...
%     'samples',Na_samples,'sampleserr',Na_sampleserr);
[Na,~,~,Na_err] = P_affine_annihilator_for_LFR(PI_a,[p_lfr;dp_lfr],'lims',pdp_lims_comp);
pcz_fhzero_report(@(pdp) Na(pdp)*PI_a(pdp), pdp_sym, Na_err, 'Annihilator Na');

%% Optimization

% Polytopes (hyperrectangles)
P_v = P_ndnorms_of_X(p_lim);
R_v = P_ndnorms_of_X(dp_lim);
PR_v = P_cartprod(P_v, R_v);

Theta_P = cellfun(@(i){ sdpvar(mx,mx,'symmetric') }, num2cell(0:np));
Q = PGenAffineMatrix(Theta_P,[1;p_sym],p_sym);
dQ = PGenAffineMatrix(Q.get_matrices,[0;dp_sym_all],pdp_sym);

% Declare other optimization variables
Lb = sdpvar(size(N.lfrtbx_obj,2), size(N.lfrtbx_obj,1), 'full');
La = sdpvar(size(Na.lfrtbx_obj,2), size(Na.lfrtbx_obj,1), 'full');
gammaSqr = sdpvar;

nr_Variables = numel([ getvariables([Q.Theta Lb]), getvariables(La + gammaSqr) ]);
pcz_dispFunction('nr. of variables: %d', nr_Variables);


CONS = [ gammaSqr >= 0 ];

for i = 1:size(P_v,1)
    p_num = P_v(i,:)';
    Nb_num = N(p_num);
    CONS = [CONS , Q(p_num) + Lb*Nb_num + Nb_num'*Lb' - eye(size(Q))*1e-5 >= 0]; %#ok<AGROW>
end

Q = Q.set_vars(pdp_sym);

for i = 1:size(PR_v,1)
    pdp_num = PR_v(i,:)';
    Na_num = Na(pdp_num);

    Q_num = Q(pdp_num);
    dQ_num = dQ(pdp_num);
    
    Qa_num = [
        He{ Ea'*Q_num*Aa } + Ea'*dQ_num*Ea + Ca'*Ca , Ea'*Q_num*Ba + Ca'*Da
        Ba'*Q_num*Ea + Da'*Ca                       , Da'*Da - gammaSqr*(Ga'*Ga)
        ];

    CONS = [CONS
        S_a'*Qa_num*S_a + La*Na_num + Na_num'*La' + eye(size(S_a,2))*1e-10 <= 0 ]; %#ok<AGROW>
    
end

Q = Q.set_vars(p_sym);

TMP_UFTXCLDbxHBtWRStETWI = pcz_dispFunctionName('Solve LMIs');

sdps = sdpsettings('solver', 'mosek', 'verbose', true, 'cachesolvers', true);
sol = optimize(CONS, gammaSqr, sdps);
pcz_feasible(sol, CONS, 'tol', 1e-6);

gamma = double(gammaSqr)^(0.5);

pcz_dispFunction
pcz_dispFunction(2,'<strong>Gamma = %g</strong> ', gamma);
pcz_dispFunction(2,'size(N)  =%4d, size(PI)   =%4d', size(N.lfrtbx_obj));
pcz_dispFunction(2,'size(Na) =%4d, size(PI_a) =%4d', size(Na.lfrtbx_obj));

Q = double(Q);
dQ = double(dQ);
Lb = double(Lb);
La = double(La);
P_cell = Q.get_matrices;
dP_cell = dQ.get_matrices;

for i = 0:np
    pcz_dispFunction_num2str(P_cell{i+1}, 'format', '%7.5g','label',sprintf('Q%d',i))    
end

pcz_dispFunction_num2str(p_lim);
pcz_dispFunction_num2str(dp_lim);
% pcz_dispFunction_num2str(Pi_indices, 'format', '%d', 'pref', ' ');
% pcz_dispFunction(msg);

pcz_dispFunctionEnd(TMP_UFTXCLDbxHBtWRStETWI);

%%
persist.stoplog;


%%


Qx = Q.set_channels([1 p_expr.'],x);

PI_x_subs = subs(PI_x(p_expr));

V = matlabFunction(x.' * PI_x_subs.' * sym(Qx) * PI_x_subs * x);

x_lim = [
    -phi_max phi_max
    -v_max v_max
    -1 1
    ];

resolution = [
    15
    15
    15
    ];

lspace = cellfun(@(o) {linspace(o{:})+eps}, num2cell(num2cell([x_lim resolution]),2));
grid = cell(1,nx);
[grid{:}] = ndgrid(lspace{:});

V(grid{:}) < 5



