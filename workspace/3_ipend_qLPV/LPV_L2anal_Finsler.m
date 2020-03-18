function LPV_L2anal_Finsler(A_fh, B_fh, C, D, K, x_lim, p_expr, p_lim, dp_lim, p_lims_comp, pdp_lims_comp, varargin)
%% LPV_L2anal_Finsler
%
%  File: LPV_L2anal_Finsler.m
%  Directory: 1_PhD_projects/22_Hinf_norm/LPV_TAC_v3_newmod
%  Author: Peter Polcz (ppolcz@gmail.com)
%
%  Created on 2018. November 21.
%

TMP_BIvQLYBqfFOwenhowBxT = pcz_dispFunctionName;

%%

args.modelname = 'model';
args.affine_P = 1;
args.skip_from_Pi = [];
args.permat = [];
args = parsepropval(args,varargin{:});

%%

O = @(varargin) zeros(varargin{:});
I = @(varargin) eye(varargin{:});
He = he;

TMP_vFNECEAeLZsYsUxvlgqL = pcz_dispFunctionName('model generation');

% Find out np, nw and nz
[nz,nx] = size(C);
[~,nw] = size(D);
np = size(p_lim,1);

% Generate symbolic variables
P_generate_symvars(nx,np,nw,nz);

A_sym = A_fh(p_cell{:});
B_sym = B_fh(p_cell{:});
Ak_sym = A_sym - B_sym*K;


%%
% We used LFR Toolbox objects (class `lfr'), which are often wrapped into a
% class `plfr' object.

% Generate symbolic lfr variables
[p_lfr,p_cell] = pcz_generateLFRStateVector('p',p_lim);
[dp_lfr,dp_cell] = pcz_generateLFRStateVector('dp',dp_lim);

p_sym = sym(plfr(p_lfr));
dp_sym = sym(plfr(dp_lfr));
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
[S_x,PI_x,iS_x,~] = P_mingen_for_LFR(PI_x,'lims',p_lims_comp);
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
[S_u,PI_u,iS_u,~] = P_mingen_for_LFR(PI_u,'lims',[-1 1]);
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
    PI_1.diff(p_cell,dp_cell)
    ];

% Lower-right block matrix of $\Pi_a$ (LFR Toolbox object)
PI_au = [ % lfr object
    PI_u.lfrtbx_obj
    PI_1 * ( F{1,2} + F{1,4}*PI_2 )
    ];

PI_a = blkdiag(PI_ax,PI_au);

% Minimal generator for PI_a
PI_a_old = plfr(PI_a,[p_lfr;dp_lfr]);
[S_a,PI_a,~,~] = P_mingen_for_LFR(PI_a_old,'lims',pdp_lims_comp);
pcz_fhzero_report(@(pdp) PI_a_old(pdp)-S_a*PI_a(pdp), pdp_sym, 'Minimal generator for PI_a(p) is OK');

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

[N,~,~,N_err] = P_affine_annihilator_for_LFR(PI_x,p_lfr,'sym',1,'lims',p_lims_comp);
pcz_fhzero_report(@(p) N(p)*PI_x(p), p_sym, N_err, 'Annihilator N');

[Na,~,~,Na_err] = P_affine_annihilator_for_LFR(PI_a,[p_lfr;dp_lfr],'lims',pdp_lims_comp);
pcz_fhzero_report(@(pdp) Na(pdp)*PI_a(pdp), pdp_sym, Na_err, 'Annihilator Na');

%% Optimization

% Polytopes (hyperrectangles)
P_v = P_ndnorms_of_X(p_lim);
R_v = P_ndnorms_of_X(dp_lim);
PR_v = P_cartprod(P_v, R_v);

Theta_P = cellfun(@(i){ sdpvar(mx,mx,'symmetric') }, num2cell(0:np));
Q = PGenAffineMatrix(Theta_P,[1;p_sym],p_sym);
dQ = PGenAffineMatrix(Q.get_matrices,[0;dp_sym],pdp_sym);

% Declare other optimization variables
Lb = sdpvar(size(N.lfrtbx_obj,2), size(N.lfrtbx_obj,1), 'full');
La = sdpvar(size(Na.lfrtbx_obj,2), size(Na.lfrtbx_obj,1), 'full');
gammaSqr = sdpvar;

nr_Variables = numel([ getvariables([Q.Theta Lb]), getvariables(La + gammaSqr) ]);
pcz_dispFunction('nr. of variables: %d', nr_Variables);


CONS = gammaSqr >= 0;

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

% sdps = sdpsettings('solver', 'mosek', 'verbose', true, 'cachesolvers', true);
sdps = sdpsettings('solver', 'SDPT3', 'verbose', true, 'cachesolvers', true);
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
    pcz_dispFunction_num2str(P_cell{i+1}, 'format', '%7.5g','name',sprintf('Q%d = ',i))    
end

pcz_dispFunction_num2str(p_lim);
pcz_dispFunction_num2str(dp_lim);
% pcz_dispFunction_num2str(Pi_indices, 'format', '%d', 'pref', ' ');
% pcz_dispFunction(msg);

pcz_dispFunctionEnd(TMP_UFTXCLDbxHBtWRStETWI);

%%

pcz_dispFunctionEnd(TMP_BIvQLYBqfFOwenhowBxT);
return

%%


Qx = Q.set_channels([1 p_expr.'],x);

PI_x_subs = subs(PI_x(p_expr));

V = matlabFunction(x.' * PI_x_subs.' * sym(Qx) * PI_x_subs * x);

resolution = [
    55
    55
    55
    ];

lspace = cellfun(@(o) {linspace(o{:})+eps}, num2cell(num2cell([x_lim resolution]),2));
xx = cell(1,nx);
[xx{:}] = ndgrid(lspace{:});

VV = V(xx{:});

pcz_2basews

alpha = 0.0031;

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
axlims = x_lim';
axis(axlims(:)),
box on
% set(gca,'LineWidth',1)

set(gca, Persist.font_axis12{:})
Persist.latexify_axis(gca,20)
Labels = Persist.latexified_labels(gca,28,'$x_1$','$x_2$','$x_3$');
Labels{1}.VerticalAlignment = 'middle';
Labels{2}.VerticalAlignment = 'bottom';

good = VV < alpha;
all = ones(size(VV));

VOLUME = prod(x_lim(:,2) - x_lim(:,1)) * ( sum(good(:))/sum(all(:)) );

display(VOLUME);

rotate3d on

end