function method4_authors_old_symbolical(modelname,A,B,C,D,p_lim,dp_lim,varargin)
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

O = @(varargin) zeros(varargin{:});
I = @(varargin) eye(varargin{:});
He = he;

TMP_vFNECEAeLZsYsUxvlgqL = pcz_dispFunctionName('model generation');

[A_sym,B_sym,C_sym,D_sym] = helper_fh2sym(A,B,C,D,p_lim);
[~,~,~,~,AC_lfr,BD_lfr,~,nx,np,nu,ny] = helper_convert(A,B,C,D,p_lim);

% 2020.04.09. (április  9, csütörtök), 09:52
% AC_lfr = minlfr(AC_lfr);
% BD_lfr = minlfr(BD_lfr);

% Generate symbolic variables
P_generate_symvars(nx,np,nu,ny);


% Find out np, nw and nz
[nz,nx] = size(C_sym);
[~,nw] = size(D_sym);
np = size(p_lim,1);

% Generate symbolic variables
P_generate_symvars(nx,np,nw,nz);

%%
% We used LFR Toolbox objects (class `lfr'), which are often wrapped into a
% class `plfr' object.

[~,AC_plfr] = P_LFR_reduction(AC_lfr,x);
[~,BD_plfr] = P_LFR_reduction(BD_lfr,w);

m1 = AC_plfr.m1;
m2 = BD_plfr.m1;
mx = nx + m1;
mu = nu + m2;

% Model matrices Fij collected in cell F{:,:}
F = cell(4,4);

[F{[1 2 3],[1 3]}] = pcz_split_matrix(AC_plfr.M,[nx,ny,m1],[nx,m1],'RowWise',0);
[F{[1 2 4],[2 4]}] = pcz_split_matrix(BD_plfr.M,[nx,ny,m2],[nu,m2],'RowWise',0);

PI_1 = sym(AC_plfr.generatePI1);
PI_x = [
    I(nx) 
    PI_1
    ];

PI_2 = sym(BD_plfr.generatePI1);
PI_u = [
    I(nu)
    PI_2
    ];

dPI_1_cell = cellfun(@(col) { jacobian(col,p)*dp }, num2cell(PI_1,1));
dPI_1 = [ dPI_1_cell{:} ];

%% Matrices of the dissipativity inequality

% Upper-left block of matrix $\Pi_a$
PI_ax = [ 
    PI_x
    PI_1 * ( F{1,1} + F{1,3}*PI_1 )
    dPI_1
    ];

% Lower-right block of matrix $\Pi_a$
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

pcz_dispFunctionEnd(TMP_vFNECEAeLZsYsUxvlgqL);

%% Optimization

% Polytopes (hyperrectangles)
P_v = P_ndnorms_of_X(p_lim);
R_v = P_ndnorms_of_X(dp_lim);
PR_v = P_cartprod(P_v, R_v);

Theta_Q = cellfun(@(i){ sdpvar(mx,mx,'symmetric') }, num2cell(0:np));
Q = PGenAffineMatrix(Theta_Q,[1;p],p);
dQ = PGenAffineMatrix(Q.get_matrices,[0;dp],pdp);

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

Q = Q.set_vars(pdp);

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

Q = Q.set_vars(p);

TMP_UFTXCLDbxHBtWRStETWI = pcz_dispFunctionName('Solve LMIs');

sdps = sdpsettings('solver', 'mosek', 'verbose', true, 'cachesolvers', true);
% sdps = sdpsettings('solver', 'SDPT3', 'verbose', true, 'cachesolvers', true);
sol = optimize(CONS, gammaSqr, sdps);
pcz_feasible(sol, CONS, 'tol', 1e-6);

gamma = double(gammaSqr)^(0.5);

Overall_Time = toc(TMP_BIvQLYBqfFOwenhowBxT);

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

pcz_dispFunction_num2str(p_lim);
pcz_dispFunction_num2str(dp_lim);
% pcz_dispFunction_num2str(Pi_indices, 'format', '%d', 'pref', ' ');
% pcz_dispFunction(msg);

store_results('Results_All.csv',modelname,0,gamma,sol.solvertime,Overall_Time,sol.info,'Polytopic a. with Finsler (old)')

pcz_dispFunctionEnd(TMP_UFTXCLDbxHBtWRStETWI);

%%

pcz_dispFunctionEnd(TMP_BIvQLYBqfFOwenhowBxT);


end