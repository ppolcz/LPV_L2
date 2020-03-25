%%
%  File: model2_ipend.m
%  Directory: 8_published/LPV_L2/workspace
%  Author: Peter Polcz (ppolcz@gmail.com)
%
%  Created on 2020. March 25. (2019b)

global SCOPE_DEPTH VERBOSE
SCOPE_DEPTH = 0;
VERBOSE = 1;

logger = Logger('results/qLPV_ipend_main-output.txt');
TMP_vcUXzzrUtfOumvfgWXDd = pcz_dispFunctionName;

RUN_ID = str2double(getenv('RUN_ID'));
if isnan(RUN_ID) || ceil(log10(RUN_ID + 1)) ~= 4
    setenv('RUN_ID', num2str(pcz_runID))
else
    setenv('RUN_ID', num2str(str2double(getenv('RUN_ID')) + 1))
end

%% Model parameters

m = 1;       % Mass of the rod [kg]
M = 2;       % Mass of the car [kg]
L = 1;       % Length of the rod [m]
g = 10;      % Gravitational acceleration [m/s^2]
k = 5;       % Friction coefficient

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
%  Compute parameter and rate bounds

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
check_stability(A_fh, B_fh, C_fh, D_fh, p_lim)

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

%%

% Imported variables to the base workspace: Q, dQ, PI_x, gamma
% method5_proposed_approach(modelname,A_fh,B_fh,C_fh,D_fh,p_lim,dp_lim,p_lims_comp,pdp_lims_comp,p_expr);


% method0_grid_author(modelname,A_fh,B_fh,C_fh,D_fh,p_lim,dp_lim,[5 5 5])
% method0_grid_author(modelname,A_fh,B_fh,C_fh,D_fh,p_lim,dp_lim,[10 10 10])
% method0_grid_author(modelname,A_fh,B_fh,C_fh,D_fh,p_lim,dp_lim,[15 15 15])
% method0_grid_author(modelname,A_fh,B_fh,C_fh,D_fh,p_lim,dp_lim,[20 20 20])

% method3_IQC_LFT_IQCToolbox(modelname,A_fh,B_fh,C_fh,D_fh,p_lim,dp_lim);
%
% Use_MinLFR = 1;
% method2_descriptor_primal(modelname,A_fh,B_fh,C_fh,D_fh,p_lim,dp_lim)
% method2_descriptor_dual(modelname,A_fh,B_fh,C_fh,D_fh,p_lim,dp_lim)
%
% opts.lpvnorm_GRID = 1;
% opts.Resolution = [5 5 5]';
% opts.lpvnorm_IQC = 1;
% opts.lpvwcgain_IQC = 1;
% LPV_L2anal_LPVTools(A_fh,B_fh,C_fh,D_fh,xw_lim,p_lim,dp_lim,opts);
%
% opts.lpvnorm_GRID = 1;
% opts.Resolution = [5 5 5]';
% opts.Resolution = [10 10 10]';
% opts.Resolution = [15 15 15]';
% opts.lpvnorm_IQC = 0;
% opts.lpvwcgain_IQC = 0;
% LPV_L2anal_LPVTools(A_fh,B_fh,C_fh,D_fh,xw_lim,p_lim,dp_lim,opts);
%
% LPV_L2anal_Finsler_old_symbolical(A_fh,B_fh,C_fh,D_fh,xw_lim,p_expr,p_lim,dp_lim,p_lims_comp,pdp_lims_comp);

%{ 
%%

TMP_oyXrvtrEzBjWNTyEyMog = pcz_dispFunctionName('Generate iso-surface of the Lyapunov/storage function');

Qx = Q.set_channels([1 p_expr.'],x);

PI_x_subs = subs(PI_x(p_expr));

V = matlabFunction(x.' * PI_x_subs.' * sym(Qx) * PI_x_subs * x);

resolution = [
    55
    55
    55
    ];

xw_lim = xw_lim(1:nx,:);
lspace = cellfun(@(o) {linspace(o{:})+eps}, num2cell(num2cell([xw_lim resolution]),2));
xx = cell(1,nx);
[xx{:}] = ndgrid(lspace{:});

VV = V(xx{:});

alpha = sdpvar;

% Hard coded (Now, I am a too lazy, to find a more adequate solution)
CONS = [
    VV([1 end],:,:) >= alpha+eps
    VV(:,[1 end],:) >= alpha+eps
    VV(:,:,[1 end]) >= alpha+eps
    ];

sol_alpha = optimize(CONS,-alpha,sdpsettings('verbose',0));
alpha = double(alpha);

pcz_dispFunction2('alpha = %g', alpha)
pcz_dispFunction2('gamma = %g', gamma)
pcz_info('|u|_L2 < M = alpha/gamma = %g', alpha/gamma)

pcz_2basews

if alpha > eps

    figure('Color', [1 1 1])
    ph = patch(isosurface(xx{:},VV,alpha));
    ph_Vertices = ph.Vertices;
    ph_Faces = ph.Faces;
    set(ph,'FaceColor','red','EdgeColor','none', 'FaceAlpha', 0.8);


    axis equal
    light('Position',[-5 0 0]),
    view([-14 15])
    xlabel $x_1$ interpreter latex
    ylabel $x_2$ interpreter latex
    zlabel $x_3$ interpreter latex

    grid on
    axlims = xw_lim';
    axis(axlims(:)),
    box on
    % set(gca,'LineWidth',1)

    set(gca, Logger.font_axis12{:})
    Logger.latexify_axis(gca,20)
    Labels = Logger.latexified_labels(gca,28,'$x_1 = v$','$x_2 = \vartheta$','$x_3 = \omega$');
    Labels{1}.VerticalAlignment = 'middle';
    Labels{2}.VerticalAlignment = 'bottom';

    good = VV < alpha;
    all = ones(size(VV));

    VOLUME = prod(xw_lim(:,2) - xw_lim(:,1)) * ( sum(good(:))/sum(all(:)) );

	pcz_dispFunction2('Approximated volume of the invariant domain: %g', VOLUME)

    rotate3d on

else

    pcz_info(0,'Iso surface completely inside X not found')

end

pcz_dispFunctionEnd(TMP_oyXrvtrEzBjWNTyEyMog);

%}

%%

pcz_dispFunctionEnd(TMP_vcUXzzrUtfOumvfgWXDd);
logger.stoplog
