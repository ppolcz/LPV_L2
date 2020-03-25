function method4_authors_old_symbolical(A_fh, B_fh, C, D, xw_lim, p_expr, p_lim, dp_lim, p_lims_comp, pdp_lims_comp, varargin)
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
args.affine_P = 0;
args.skip_from_Pi = [];
args.permat = [];
%%
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


%%
% We used LFR Toolbox objects (class `lfr'), which are often wrapped into a
% class `plfr' object.

% Generate symbolic lfr variables
[p_lfr,p_cell] = pcz_generateLFRStateVector('p',p_lim);
[dp_lfr,dp_cell] = pcz_generateLFRStateVector('dp',dp_lim);

p_sym = sym(plfr(p_lfr));
dp_sym = sym(plfr(dp_lfr));
pdp_sym = [ p_sym ; dp_sym ];


% LFR realization of matrix A(p):
A_lfr = minlfr(A_fh(p_cell{:}));

% LFR realization of matrix B(p):
B_lfr = minlfr(B_fh(p_cell{:}));

A_sol = P_LFR_reduction(A_lfr,x);
A_plfr = A_sol.LFR_reduced;

B_sol = P_LFR_reduction(B_lfr,w);
B_plfr = B_sol.LFR_reduced;

% -----------------------

% Model matrices Fij collected in cell F{:,:}
F = cell(4,4);

%%%
% LFR of matrix A(p) = (F11 F12)*PI_x, where PI_x = [ eye(3) ; PI_1 ]
% 
% Subscript |_x| is not used in the manuscript. |PI_x|, |S_x| and |iS_x| in
% the script corresponds to $\Pi$, $S$ and $S^{-1}$ in the manuscript.
[F{1,1},F{1,3},F{3,1},F{3,3},PI_x,PI_1] = deal(A_plfr.A, A_plfr.B, A_plfr.C, A_plfr.D, A_plfr.PI_b, A_plfr.PI);
pcz_symzero_report(A_sym - [ F{1,1} F{1,3} ] * PI_x,sprintf('Generator form (m = %d) of A(p)-B(p)*K is OK', size(PI_1,1)));
m1 = size(F{1,3},2);
mx = nx + m1;

% C(x,p)
F{2,1} = C;
F{2,3} = O(nz,m1);

% -----------------------

% B(x,p)
[F{1,2},F{1,4},F{4,2},F{4,4},PI_u,PI_2] = deal(B_plfr.A, B_plfr.B, B_plfr.C, B_plfr.D, B_plfr.PI_b,B_plfr.PI);
pcz_symzero_report(B_sym - [ F{1,2} F{1,4} ] * PI_u,'Generator form of B(p) is OK');
m2 = size(F{1,4},2);
mu = nw + m2;

% D(x,p)
F{2,2} = D;
F{2,4} = O(nz,m2);

%% Matrices of the dissipativity inequality

% Upper-left block matrix of $\Pi_a$
PI_ax = [ 
    PI_x
    PI_1 * ( F{1,1} + F{1,3}*PI_1 )
    diff(PI_1,p1)*dp1 + diff(PI_1,p2)*dp2 + diff(PI_1,p3)*dp3
    ];

% Lower-right block matrix of $\Pi_a$
PI_au = [
    PI_u
    PI_1 * ( F{1,2} + F{1,4}*PI_2 )
    ];

PI_a = blkdiag(PI_ax,PI_au);

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

N = P_affine_annihilator(PI_x*x,xp,p);
N = N.set_subsvars(p);

Na = P_affine_annihilator(PI_a*xw,xwpdp,pdp);
Na = Na.set_subsvars(pdp);

%% Optimization

% Polytopes (hyperrectangles)
P_v = P_ndnorms_of_X(p_lim);
R_v = P_ndnorms_of_X(dp_lim);
PR_v = P_cartprod(P_v, R_v);

Theta_Q = cellfun(@(i){ sdpvar(mx,mx,'symmetric') }, num2cell(0:np));
Q = PGenAffineMatrix(Theta_Q,[1;p_sym],p_sym);
dQ = PGenAffineMatrix(Q.get_matrices,[0;dp_sym],pdp_sym);

% Declare other optimization variables
Lb = sdpvar(size(N,2), size(N,1), 'full');
La = sdpvar(size(Na,2), size(Na,1), 'full');
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
        Qa_num + La*Na_num + Na_num'*La' + eye(size(Qa_num,2))*1e-10 <= 0 ]; %#ok<AGROW>
    
end

Q = Q.set_vars(p_sym);

TMP_UFTXCLDbxHBtWRStETWI = pcz_dispFunctionName('Solve LMIs');

sdps = sdpsettings('solver', 'mosek', 'verbose', true, 'cachesolvers', true);
% sdps = sdpsettings('solver', 'SDPT3', 'verbose', true, 'cachesolvers', true);
sol = optimize(CONS, gammaSqr, sdps);
pcz_feasible(sol, CONS, 'tol', 1e-6);

gamma = double(gammaSqr)^(0.5);

pcz_dispFunction
pcz_dispFunction(2,'<strong>Gamma = %g</strong> ', gamma);
pcz_dispFunction(2,'size(N)  =%4d, size(PI)   =%4d', size(N));
pcz_dispFunction(2,'size(Na) =%4d, size(PI_a) =%4d', size(Na));

Q = double(Q);
dQ = double(dQ);
Lb = double(Lb);
La = double(La);
Q_cell = Q.get_matrices;
dP_cell = dQ.get_matrices;

pcz_2basews(Q,dQ,PI_x,gamma)

for i = 0:np
    pcz_dispFunction_num2str(Q_cell{i+1}, 'format', '%7.5g','name',sprintf('Q%d = ',i))    
end

pcz_dispFunction_num2str(xw_lim);
pcz_dispFunction_num2str(p_lim);
pcz_dispFunction_num2str(dp_lim);
% pcz_dispFunction_num2str(Pi_indices, 'format', '%d', 'pref', ' ');
% pcz_dispFunction(msg);

store_results('Results_All.csv',xw_lim,[],gamma,sol.solvertime,sol.info,'Polytopic a. with Finsler (old)',0)

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

lspace = cellfun(@(o) {linspace(o{:})+eps}, num2cell(num2cell([xw_lim resolution]),2));
xx = cell(1,nx);
[xx{:}] = ndgrid(lspace{:});

VV = V(xx{:});

% pcz_2basews

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
axlims = xw_lim';
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

VOLUME = prod(xw_lim(:,2) - xw_lim(:,1)) * ( sum(good(:))/sum(all(:)) );

display(VOLUME);

rotate3d on

end