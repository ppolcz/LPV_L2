%% 
%  File: model4_randLFR.m
%  Directory: 8_published/LPV_L2/workspace
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2020. March 28. (2019b)
%

%%
% Automatically generated stuff

global SCOPE_DEPTH VERBOSE 
SCOPE_DEPTH = 0;
VERBOSE = 1;

try c = evalin('caller','persist'); catch; c = []; end
persist = Persist(mfilename('fullpath'), c); clear c; 
persist.backup();
%clear persist

%%

O = @(varargin) zeros(varargin{:});
I = @(varargin) eye(varargin{:});
He = he;

nx = 2;
np = 1;
nu = 1;
ny = 1;

r = [0 1];
m1 = sum(r);
m = nx + m1;

F = cell(3,2);
F_M = [
    -7  1 , -4
     1  1 ,  2
    ..........
     3 -2 ,  3
    ..........
    -2  0 ,  0
    ];
[F{:}] = pcz_split_matrix(F_M,[nx nu m1],[nx m1],'RowWise',0);

C = [1 1];
D = 0;

p_lim = repmat([-1,1],[np,1]);
dp_lim = repmat([-1,1],[np,1]);

pdp_lim = [ p_lim ; dp_lim ];

P_generate_symvars(nx,np,nu,ny);
[p_lfr,p_lfr_cell] = pcz_generateLFRStateVector('p',p_lim);
[dp_lfr,dp_lfr_cell] = pcz_generateLFRStateVector('dp',p_lim);

Delta = cellfun(@(ri,pi) { I(ri)*pi }, num2cell(r), [1 p_lfr_cell]);
Delta = blkdiag(Delta{:});

PI_1 = plfr(inv( I(m1) - Delta*F{3,2} ) * Delta*F{3,1});
PI = plfr([ I(nx) ; PI_1 ]);

A_lfr = PI' * [ F{1,1} F{1,2} ]';
B_lfr = PI' * [ F{2,1} F{2,2} ]';

%{

% Explicit symbolical value:
Delta_sym = sym(plfr(Delta));
PI_1_sym = ( I(m1) - Delta_sym*F{3,2} ) \ Delta_sym*F{3,1};
PI_sym = [ I(nx) ; PI_1_sym ];

A_sym = ( [F{1,1} F{1,2}] * PI_sym )';
B_sym = ( [F{2,1} F{2,2}] * PI_sym )';

%}

PI_d = [
    PI
    PI_1*F{1,1}
    PI_1*F{1,2}*PI_1
    diff(PI_1,p_lfr_cell,dp_lfr_cell,pdp)
    ];

md = size(PI_d,1);

Ad = [
    F{1,1}   F{1,2}   O(nx,m1) O(nx,m1) O(nx,m1)
    O(m1,nx) O(m1,m1) I(m1)    I(m1)   -I(m1)
    ];

Ed = [
    I(nx)    O(nx,m1) O(nx,m1) O(nx,m1) O(nx,m1)
    O(m1,nx) I(m1)    O(m1,m1) O(m1,m1) O(m1,m1)
    ];

Bd = [
    F{2,1}   F{2,2}   O(nu,m1) O(nu,m1) O(nu,m1)
    ];

%% Construct annihilators

[N,~,~,N_err] = P_affine_annihilator_for_LFR(PI,p_lfr,'lims',p_lim);
pcz_fhzero_report(@(p) N(p)*PI(p), p, N_err*10, 'Annihilator N');

[Nd,~,~,Nd_err] = P_affine_annihilator_for_LFR(PI_d,[p_lfr;dp_lfr],'lims',pdp_lim);
pcz_fhzero_report(@(pdp) Nd(pdp)*PI_d(pdp), pdp, Nd_err, 'Annihilator Nd');

%%

% Polytopes (hyperrectangles)
P_v = P_ndnorms_of_X(p_lim);
R_v = P_ndnorms_of_X(dp_lim);
PR_v = P_cartprod(P_v, R_v);

nP = size(P_v,1);
nR = size(R_v,1);
nPR = size(PR_v,1);

ind = 1+0*(0:np);

Theta_Q = sdpvar(m*ind);
Theta_H = sdpvar(nu*ind,m*ind,'full');

Q = PAffineMatrix(Theta_Q,[1;p],p);
dQ = Q.set_channels([0;dp],pdp);
H = PAffineMatrix(Theta_H,[1;p],pdp);

% Declare other optimization variables
L = sdpvar(size(N,2), size(N,1), 'full');
Ld = sdpvar(size(Nd,2), size(Nd,1), 'full');

CONS = [];

for i = 1:nP
    p_num = P_v(i,:)';
    
    LMI_1 = Q(p_num) + He{ L*N(p_num) } - blkdiag(I(nx),O(m-nx))*1e-5;
    
    CONS = [CONS , LMI_1 >= 0]; %#ok<AGROW>
end

Qd = He{ Ed'*Q*Ad - Ed'*H'*Bd } - Ed'*dQ*Ed;
Qd = Qd.set_subsvars(pdp);

for i = 1:nPR
    pdp_num = PR_v(i,:)';
    
    LMI_2 = Qd(pdp_num) + He{ Ld*Nd(pdp_num) } + blkdiag(I(nx),O(md-nx))*1e-5;

    CONS = [CONS , LMI_2 <= 0]; %#ok<AGROW>
end

TMP_UFTXCLDbxHBtWRStETWI = pcz_dispFunctionName('Solve LMIs');

sdps = sdpsettings('solver', 'lmilab', 'verbose', true, 'cachesolvers', true);
sol = optimize(CONS, [], sdps);
pcz_feasible(sol, CONS, 'tol', 1e-6);
pcz_feasible_Finsler(Q,L,N,p_lim,'title','Dissipativity relation')
pcz_feasible_Finsler(-Qd,-Ld,Nd,pdp_lim,'title','Dissipativity relation')

Q = double(Q);
L = double(L);
Qd = double(Qd);
Ld = double(Ld);

LMI_2 = He{ Ld*Nd } + plfr(Qd,pdp);
helper_check_stability(LMI_2,pdp_lim)

LMI_2 = PI_d' * plfr(double(Qd)) * PI_d;
helper_check_stability(LMI_2,pdp_lim)

%% Compute state feedback gain

QQ = PI' * plfr(double(Q)) * PI;
HH = plfr(double(H)) * PI;

% H = K*Q  ==>  K = H/Q
KK = HH / QQ;

Ak_lfr = A_lfr - B_lfr*KK;
helper_check_stability(Ak_lfr,p_lim)

%%

A = Ak_lfr;
B = B_lfr;

bases = [
    1
    p1
    p1^2
    p1^3
    p1^4
    p1^5
    ];

bases_Jac = jacobian(bases,p);

modelname = 'simple2D_1x1_LFR';

method1_RCT(modelname,A,B,C,D,p_lim)

method0_grid_LPVTools(modelname,A,B,C,D,p_lim,dp_lim,bases,bases_Jac,[5 5 5]);
% 
% Greedy grid 
method0_grid_Wu1995(modelname,A,B,C,D,p_lim,dp_lim,bases,bases_Jac,'res_max',5');
% method0_grid_Wu1995(modelname,A,B,C,D,p_lim,dp_lim,bases,bases_Jac,'res_max',15');
% 
% % As proposed by Wu (1995,1996)
% method0_grid_Wu1995(modelname,A,B,C,D,p_lim,dp_lim,bases,bases_Jac,'res_max',1000,'delta',1e-4,'T',10000);
% method0_grid_Wu1995(modelname,A,B,C,D,p_lim,dp_lim,bases,bases_Jac,'res_max',1000,'delta',1e-3,'T',1000);
% method0_grid_Wu1995(modelname,A,B,C,D,p_lim,dp_lim,bases,bases_Jac,'res_max',1000,'delta',1e-2,'T',100);
% method0_grid_Wu1995(modelname,A,B,C,D,p_lim,dp_lim,bases,bases_Jac,'res_max',1000,'delta',1e-1,'T',100);
% method0_grid_Wu1995(modelname,A,B,C,D,p_lim,dp_lim,bases,bases_Jac,'res_max',1000,'delta',1e-2,'T',10);
% method0_grid_Wu1995(modelname,A,B,C,D,p_lim,dp_lim,bases,bases_Jac,'res_max',1000,'delta',1e-1,'T',10);
% 
method2_descriptor_primal(modelname,A,B,C,D,p_lim,dp_lim)
method2_descriptor_dual(modelname,A,B,C,D,p_lim,dp_lim,1)
method2_descriptor_dual(modelname,A,B,C,D,p_lim,dp_lim,0)
% 
% IQC/LFT approaches for LPV with rate-bounded parameters
method3_IQC_LFT_IQCToolbox(modelname,A,B,C,D,p_lim,dp_lim);
method3_IQC_LFT_LPVTools(modelname,A,B,C,D,p_lim,dp_lim);

method4_authors_old_symbolical(modelname,A,B,C,D,p_lim,dp_lim);

% Imported variables to the base workspace: Q, dQ, PI_x, gamma
method5_proposed_approach(modelname,A,B,C,D,p_lim,dp_lim,p_lim,pdp_lim,[]);

%%
persist.stoplog;
