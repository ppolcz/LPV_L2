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

rnd = @(varargin) randn(varargin{:});

nx = 3;
np = 2;
nu = 2;
ny = 1;

r = [0 2 1];
m1 = sum(r);
m = nx + m1;

p_lim = repmat([-1,1],[np,1]);
dp_lim = repmat([-0.1,0.1],[np,1]);

pdp_lim = [ p_lim ; dp_lim ];

P_generate_symvars(nx,np,nu,ny);
[p_lfr,p_lfr_cell] = pcz_generateLFRStateVector('p',p_lim,dp_lim);
[dp_lfr,dp_lfr_cell] = pcz_generateLFRStateVector('dp',p_lim);
pdp_lfr_cell = [ p_lfr_cell dp_lfr_cell ];

F = cell(3,2);

[F{:}] = pcz_split_matrix(rnd(nx+nu+m1,nx+m1),[nx nu m1],[nx m1],'RowWise',0);
C = rnd(ny,nx);
D = rnd(ny,nu);

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

% %% Construct annihilators

[N,~,~,N_err] = P_affine_annihilator_for_LFR(PI,p_lfr,'lims',p_lim);
pcz_fhzero_report(@(p) N(p)*PI(p), p, N_err*10, 'Annihilator N');

[Nd,~,~,Nd_err] = P_affine_annihilator_for_LFR(PI_d,[p_lfr;dp_lfr],'lims',pdp_lim);
pcz_fhzero_report(@(pdp) Nd(pdp)*PI_d(pdp), pdp, Nd_err, 'Annihilator Nd');

% %%

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

    LMI_1 = Q(p_num) + He{ L*N(p_num) } - I(m)*1e-5;

    CONS = [CONS , LMI_1 >= 0]; %#ok<AGROW>
end

Qd = He{ Ed'*Q*Ad - Ed'*H'*Bd } - Ed'*dQ*Ed;
Qd = Qd.set_subsvars(pdp);

for i = 1:nPR
    pdp_num = PR_v(i,:)';

    LMI_2 = Qd(pdp_num) + He{ Ld*Nd(pdp_num) } + I(md)*1e-5;

    CONS = [CONS , LMI_2 <= 0]; %#ok<AGROW>
end

TMP_UFTXCLDbxHBtWRStETWI = pcz_dispFunctionName('Solve LMIs');

sdps = sdpsettings('solver', 'mosek', 'verbose', true, 'cachesolvers', true);
sol = optimize(CONS, [], sdps);
pcz_feasible(sol, CONS, 'tol', 1e-6);
pcz_feasible_Finsler(Q,L,N,p_lim,'title','Dissipativity relation')
pcz_feasible_Finsler(-Qd,-Ld,Nd,pdp_lim,'title','Dissipativity relation')

H = double(H);
Q = double(Q);
L = double(L);
Qd = double(Qd);
Ld = double(Ld);

H_lfr = plfr(H,p_lfr_cell);
Q_lfr = plfr(Q,p_lfr_cell);
Qd_lfr = plfr(Qd,pdp_lfr_cell);

%%

LMI_2 = He{ Ld*Nd } + Qd_lfr;
helper_check_stability(LMI_2,pdp_lim)

LMI_2 = PI_d' * Qd_lfr * PI_d;
helper_check_stability(LMI_2,pdp_lim)

%%

QQ = PI' * Q_lfr * PI;
HH = H_lfr * PI;

% H = K*Q  ==>  K = H/Q
KK = HH / QQ;

%%

Ak_lfr = A_lfr - B_lfr*KK;


helper_check_stability(Ak_lfr,p_lim)

% return

%%

A = Ak_lfr;
B = B_lfr;

bases = [
    1
    p1
    p1^2
    p2
    p2^2
    ];

bases_Jac = jacobian(bases,p);

modelname = sprintf('random_%dD_%dx%d_LFR',nx,ny,nu);

method1_RCT(modelname,A,B,C,D,p_lim)

% method0_grid_LPVTools(modelname,A,B,C,D,p_lim,dp_lim,bases,bases_Jac,[5 5 5]);
% 
% Greedy grid 
% method0_grid_Wu1995(modelname,A,B,C,D,p_lim,dp_lim,bases,bases_Jac,'res_max',5');
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
% method2_descriptor_primal(modelname,A,B,C,D,p_lim,dp_lim)
% method2_descriptor_dual(modelname,A,B,C,D,p_lim,dp_lim,1)
% method2_descriptor_dual(modelname,A,B,C,D,p_lim,dp_lim,0)
% 
% % IQC/LFT approaches for LPV with rate-bounded parameters
% method3_IQC_LFT_IQCToolbox(modelname,A,B,C,D,p_lim,dp_lim);
% method3_IQC_LFT_LPVTools(modelname,A,B,C,D,p_lim,dp_lim);
% 
% method4_authors_old_symbolical(modelname,A,B,C,D,p_lim,dp_lim);

% Imported variables to the base workspace: Q, dQ, PI_x, gamma
method5_proposed_approach(modelname,A,B,C,D,p_lim,dp_lim,p_lim,pdp_lim,[]);

%%
persist.stoplog;




%%
%{

F{1,1} = [
     6.547007489022763e-01    -2.240119862734167e+00    -3.479254962059620e-01     8.944613154736245e-02    -5.603617495969332e-01
     2.868112010418336e-01    -2.356442799588965e-01     5.073296875899195e-01     2.782586115601847e+00     4.334415720667913e-01
     1.418357597818456e+00    -1.253171284466545e+00    -1.222167502701937e+00     1.813291347063874e-01     4.707563518173649e-01
    -1.118931815980435e+00    -1.030196066545677e+00     2.156430125953823e-02    -5.969451513596548e-01    -8.425704527608475e-01
    -7.468767781339673e-02     5.928356915286037e-01    -8.857083344913315e-01     1.002224384479816e+00     3.759190094272149e-01
    ];

F{2,1} = [
     1.084119485612190e+00     1.090252577340786e+00    -1.916228167112975e-01     9.248121003585645e-01    -1.121856471186102e+00
    -1.454956239027367e-02    -2.094614038416507e-01     6.895145533682722e-01    -1.521679472588687e+00     8.248895905675911e-01
     6.429255414390348e-01    -5.192212120606393e-01     1.959331995574537e+00    -6.400853377089221e-01    -2.138467620286843e-01
    ];

F{3,1} = [
    -1.619869242827445e+00     1.802246084853087e+00    -4.520414245013743e-01    -4.480152096763123e-01     7.200252243939719e-01
     3.487475292046774e-01    -8.612515661978518e-02     5.332546267796872e-01    -2.536666831871017e-01    -9.934682510089320e-01
     1.083413492479740e+00     3.806104303132576e-01    -4.570554443207401e-01     1.049093173475492e-01     7.818648172900673e-01
    ];

F{1,2} = [
    -2.762671926502525e-01     9.420729181773211e-01    -9.742493462098838e-01
    -7.264128972024417e-02    -8.749846035850006e-02    -9.795011891645293e-01
    -2.717237474273746e-01     1.192776006286231e+00    -9.214520844692828e-01
     6.712409704845821e-01    -4.981682153209956e-01     8.231383883151274e-01
     1.802719740002845e+00     4.375462697296674e-01     6.565076507505641e-03
    ];

F{2,2} = [
     8.960985173123124e-01     3.300441888842542e-01     3.741327817037007e-01
     2.944991646799471e-02    -1.400141860139974e+00    -5.551952188032477e-01
    -7.475762176852034e-01    -5.218810840534620e-01    -1.246384660283821e-01
    ];

F{3,2} = [
     8.404722698175919e-02     1.004918393567770e+00     8.369742504322356e-01
    -9.976195945410180e-01    -1.137789385757429e+00    -1.287476203065973e+00
    -8.735719145174653e-02    -1.532479975942847e-01    -7.585127542204142e-01
    ];
%}