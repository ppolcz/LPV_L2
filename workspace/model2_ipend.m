%%
%  File: model2_ipend.m
%  Directory: 8_published/LPV_L2/workspace
%  Author: Peter Polcz (ppolcz@gmail.com)
%
%  Created on 2020. March 25. (2019b)

G_reset
P_init(11)

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

%% Model parameters

m = 1;       % Mass of the rod [kg]
M = 2;       % Mass of the car [kg]
L = 1;       % Length of the rod [m]
g = 10;      % Gravitational acceleration [m/s^2]
k = 5;       % Friction coefficient [kg/s]

% Moment of inertia of a rod of length 2L and mass m, rotating about one
% end. This expression assumes that the rod is an infinitely thin (but
% rigid) wire.
I = 4*m*L^2/3;

pcz_dispFunction2('Model parameters:')
pcz_dispFunction_scalar(m,M,L,g,k,I)

%%%
% Generate x1,x2,x3; p1,p2,p3; w1; z1
P_generate_symvars(3,3,1,1);
veloc = x1;
theta = x2;
omega = x3;

p_expr = [
    sin(theta)/theta  % = p1
    cos(theta)        % = p2
    omega*sin(theta)  % = p3
    ];
p_expr_cell = num2cell(p_expr);

%% Model matrices

sigma1 = @(p2) (I+m*L^2)*(m+M) - m^2*L^2*p2^2;

A_fh = @(p1,p2,p3) [
    -k*(I+m*L^2) , -m^2*L^2*g*p1*p2 , -m*L*(I+m*L^2)*p3
     0           ,  0               ,  sigma1(p2)
    -k*m*L*p2    , -(m+M)*m*g*L*p1  , -m^2*L^2*p2*p3
    ] / sigma1(p2);

B_fh = @(p1,p2,p3) [
    I+m*L^2
    0
    m*L*p2
    ] / sigma1(p2);

C_fh = [ 0 1 0 ];

D_fh = 0;

% Needed to compute the rate bounds of p3
Ax_sym = A_fh(p_expr_cell{:});
Bx_sym = B_fh(p_expr_cell{:});

%% Bound the state and input space
%  Compute the parameter and rate bounds

% Bound the state variables and input
veloc_max = 2;   % meters/seconds
theta_max = 0.4; % radians
omega_max = 1.2; % radians/seconds
w_max = 1;       % bound the disturbance input

modelname = sprintf('ipend(%g,%g,%g,%g,%g;%g,%g,%g;%g)', m,M,L,g,k, ...
    veloc_max, theta_max, omega_max, w_max);

pcz_dispFunction2('State and input bounds:')
pcz_dispFunction_scalar(veloc_max, theta_max, omega_max, w_max);

xw_lim = [
    -veloc_max  veloc_max
    -theta_max  theta_max
    -omega_max  omega_max
    -w_max      w_max
    ];

% Compute numerically the bounds for dp3
resolution = [
    45
    45
    45
    2
    ];
lspace = cellfun(@(o) {linspace(o{:})+eps}, num2cell(num2cell([xw_lim resolution]),2));
grid = cell(1,size(xw_lim,1));
[grid{:}] = ndgrid(lspace{:});

dx3_fh = matlabFunction(expand(simplify(Ax_sym(3,:)*x + Bx_sym(3,:)*w)),'vars',[x;w]);
dx3_grid = dx3_fh(grid{:});

p1_lim = [ sin(theta_max) / theta_max 1 ];
p2_lim = [ cos(theta_max) 1 ];
p3_lim = omega_max * [ -sin(theta_max) sin(theta_max) ];

dp1_fn = @(x1,x3) ( x1*cos(x1) - sin(x1) ) / x1^2 * x3;
dp1_lim = [ dp1_fn(theta_max,omega_max) dp1_fn(-theta_max,omega_max) ];
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
helper_check_stability(A_fh, B_fh, C_fh, D_fh, p_lim)

% Give some nice domains, where the rational matrices are well-defined and
% on which the numerical computations will be performed (minimal generator
% and annihilator computation).
p_lims_comp = [
    -1 1
    -1 1
    -1 1
    ];
pdp_lims_comp = [
    p_lims_comp
    -1 1
    -1 1
    -1 1
    ];

%% Basis functions for the grid-based methods

bases = [
                   1
                  p1
                  p2
                  p3
 -(p1*p2)/(p2^2 - 7)
      -p1/(p2^2 - 7)
      -p2/(p2^2 - 7)
     p2^2/(p2^2 - 7)
 -(p2*p3)/(p2^2 - 7)
      -p3/(p2^2 - 7)
    ];

bases_Jac = jacobian(bases,p);


%%

method0_grid_ltiwc_Hinf(modelname,A_fh,B_fh,C_fh,D_fh,p_lim)

% method1_RCT(modelname,A_fh,B_fh,C_fh,D_fh,p_lim)

% method0_grid_LPVTools(modelname,A_fh,B_fh,C_fh,D_fh,p_lim,dp_lim,bases,bases_Jac,[5 5 5]);
% 
% % Greedy grid 
% method0_grid_Wu1995(modelname,A_fh,B_fh,C_fh,D_fh,p_lim,dp_lim,bases,bases_Jac,'res_max',5);
% method0_grid_Wu1995(modelname,A_fh,B_fh,C_fh,D_fh,p_lim,dp_lim,bases,bases_Jac,'res_max',15);
% 
% % As proposed by Wu (1995,1996)
% method0_grid_Wu1995(modelname,A_fh,B_fh,C_fh,D_fh,p_lim,dp_lim,bases,bases_Jac,'res_max',5,'delta',1e-6,'T',1e6);
% method0_grid_Wu1995(modelname,A_fh,B_fh,C_fh,D_fh,p_lim,dp_lim,bases,bases_Jac,'res_max',5,'delta',1e-5,'T',100000);
% method0_grid_Wu1995(modelname,A_fh,B_fh,C_fh,D_fh,p_lim,dp_lim,bases,bases_Jac,'res_max',5,'delta',1e-4,'T',10000);
% method0_grid_Wu1995(modelname,A_fh,B_fh,C_fh,D_fh,p_lim,dp_lim,bases,bases_Jac,'res_max',5,'delta',1e-3,'T',1000);
% method0_grid_Wu1995(modelname,A_fh,B_fh,C_fh,D_fh,p_lim,dp_lim,bases,bases_Jac,'res_max',5,'delta',1e-2,'T',100);
% method0_grid_Wu1995(modelname,A_fh,B_fh,C_fh,D_fh,p_lim,dp_lim,bases,bases_Jac,'res_max',5,'delta',1e-1,'T',100);
% method0_grid_Wu1995(modelname,A_fh,B_fh,C_fh,D_fh,p_lim,dp_lim,bases,bases_Jac,'res_max',5,'delta',1e-2,'T',10);
% method0_grid_Wu1995(modelname,A_fh,B_fh,C_fh,D_fh,p_lim,dp_lim,bases,bases_Jac,'res_max',5,'delta',1e-1,'T',10);
% 
% method2_descriptor_primal(modelname,A_fh,B_fh,C_fh,D_fh,p_lim,dp_lim)
% method2_descriptor_dual(modelname,A_fh,B_fh,C_fh,D_fh,p_lim,dp_lim,1)
% method2_descriptor_dual(modelname,A_fh,B_fh,C_fh,D_fh,p_lim,dp_lim,0)
% 
% % IQC/LFT approaches for LPV with rate-bounded parameters
% method3_IQC_LFT_IQCToolbox(modelname,A_fh,B_fh,C_fh,D_fh,p_lim,dp_lim);fi
% method3_IQC_LFT_LPVTools(modelname,A_fh,B_fh,C_fh,D_fh,p_lim,dp_lim);
% 
% method4_authors_old_symbolical(modelname,A_fh,B_fh,C_fh,D_fh,p_lim,dp_lim);
% 
% % Imported variables to the base workspace: Q, dQ, PI_x, gamma
% method5_proposed_approach(modelname,A_fh,B_fh,C_fh,D_fh,p_lim,dp_lim,p_lims_comp,pdp_lims_comp,p_expr);
% helper_generate_isosurface_method5(Q,PI_x,gamma,p_expr,xw_lim(1:nx,:),p_lim,dp_lim,[1,1,1]*101);

%% 
%{

global LPVTools_vars

IQCinfoS = LPVTools_vars.IQCinfoS

omega = LPVTools_vars.omega;
[~,index] = min(abs(omega - 2*pi));
omega = omega(index);

Psi = cell(1,np);
[Psi{:}] = LPVTools_vars.IQCinfoS.PsiPi;

Ts_sym = sym([])
Ts = [];
Us = [];
Vs = [];
syms s
for i = 1:np
    pl = LPVTools_vars.polelist.(['p' num2str(i)]);
    pl(isinf(pl)) = [];
    pl = abs(1./(1i*omega-pl)');
    % pl = reshape(pl, [1,1,numel(pl)]);
        
    z = sym('z',size(pl));
    
    Psi_i1 = abs(Psi{i}{1}(:,:,index));
    Psi_i2 = abs(Psi{i}{2}(:,:,index));
    
    Psi_i = [
        Psi_i1(1:IQCinfoS(i).ExtBlkDim(1),:)
        Psi_i2(1:numel(pl),:)
        Psi_i1(IQCinfoS(i).ExtBlkDim(1)+1:end,:)
        Psi_i2(numel(pl)+1:end,:)
    ];
    
    Tsi = Psi_i1(1:IQCinfoS(i).ExtBlkDim(1),1:IQCinfoS(i).OrigBlkDim(1));
    Tsi_sym = sym(round(Tsi,10));
    Tsi_sym(abs(Tsi - 1) > 1e-10) = 0;
    for j = 1:numel(pl)
        Tsi_sym(abs(Tsi - pl(j)) < 1e-10) = 1/(s - z(j));
    end

    Ts_sym = blkdiag(Ts_sym, Tsi_sym);
    Ts = blkdiag(Ts, Tsi);
    Us = blkdiag(Us, Psi_i2(1:numel(pl),1:IQCinfoS(i).OrigBlkDim(1)));    
    Vs = blkdiag(Vs, Psi_i1(IQCinfoS(i).ExtBlkDim(1)+1:end,2*IQCinfoS(i).OrigBlkDim(2)+1:end));
end

%}
%%

pcz_dispFunctionEnd(TMP_vcUXzzrUtfOumvfgWXDd);
logger.stoplog

return

%%

% Generate grid
lspace = cellfun(@(o) {linspace(o{:})}, num2cell(num2cell([p_lim res_used]),2));
pp = cell(1,np);
[pp{:}] = ndgrid(lspace{:});
pp = cellfun(@(a) {a(:)'}, pp);
pp = vertcat(pp{:});


norm(ss(A_fh(pwc{:}),B_fh(pwc{:}),C_fh,D_fh),Inf)

