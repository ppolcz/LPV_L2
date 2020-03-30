function method5_proposed_approach(modelname, A_fh, B_fh, C_fh, D_fh, p_lim, dp_lim, p_lims_comp, pdp_lims_comp, varargin)
%% 
%  File: method5_proposed_approach.m
%  Directory: 8_published/LPV_L2/workspace
%  Author: Peter Polcz (ppolcz@gmail.com) 
% 
%  Created on 2018. November 21.
%  Major review on 2020. March 25. (2019b)
%

TMP_BIvQLYBqfFOwenhowBxT = pcz_dispFunctionName;

%%

O = @(varargin) zeros(varargin{:});
I = @(varargin) eye(varargin{:});
He = he;

TMP_vFNECEAeLZsYsUxvlgqL = pcz_dispFunctionName('model generation');

[~,~,~,~,AC_lfr,BD_lfr] = helper_convert(A_fh, B_fh, C_fh, D_fh, p_lim);

% Dimensions
nx = size(AC_lfr,2);
ny = size(AC_lfr,1) - nx;
nu = size(BD_lfr,2);
np = size(p_lim,1);

% Generate symbolic variables
P_generate_symvars(nx,np,nu,ny);

% Generate lfr variables
[p_lfr,p_cell] = pcz_generateLFRStateVector('p',p_lim);
[dp_lfr,dp_cell] = pcz_generateLFRStateVector('dp',dp_lim);


%%
% We used LFR Toolbox objects (class `lfr'), which are often wrapped into a
% class `plfr' object.

% Model matrices Fij collected in cell F{:,:}
F = cell(4,4);

% -----------------------

%%%
%  LFR of matrix [ A(p) ] = [ F11 F13 ] * PI_x, where PI_x = [ eye(nx) ; PI_1 ]
%                [ C(p) ]   [ F21 F23 ]
% 
AC_plfr = plfr(AC_lfr);
[F{1,1},F{1,3},F{3,1},F{3,3},PI_x,PI_1] = deal(AC_plfr.A, AC_plfr.B, AC_plfr.C, AC_plfr.D, AC_plfr.generatePI, AC_plfr.generatePI1);
pcz_lfrzero_report(AC_lfr - [ F{1,1} F{1,3} ] * PI_x, sprintf('Generator form (m = %d) of [ A(p) ; C(p) ] is OK', PI_1.ny));

% Minimal generator for PI_x
[S_x,PI_x,iS_x,~] = P_mingen_for_LFR(PI_x.set_vars(p),'lims',p_lims_comp);
[F{1,1},F{1,3},F{3,1},F{3,3}] = pcz_split_matrix([F{1,1} F{1,3} ; F{3,1} F{3,3}] * S_x, [nx+nu Inf], [nx Inf]);
pcz_lfrzero_report(AC_lfr - [ F{1,1} F{1,3} ] * PI_x,'Minimal generator form of [ A(p) ; C(p) ] is OK');

% Reload PI_1 corresponding to the minimal generator PI_x.
PI_1 = plfr(iS_x(nx+1:end,nx+1:end) * PI_1); % plfr object

% -----------------------

%%%
%  LFR of matrix [ B(p) ] = [ F12 F14 ] * PI_u, where PI_u = [ eye(nu) ; PI_2 ]
%                [ D(p) ]   [ F22 F24 ]
% 
BD_plfr = plfr(BD_lfr);
[F{1,2},F{1,4},F{4,2},F{4,4},PI_u,PI_2] = deal(BD_plfr.A, BD_plfr.B, BD_plfr.C, BD_plfr.D, BD_plfr.generatePI,BD_plfr.generatePI1);
pcz_lfrzero_report(BD_lfr - [ F{1,2} F{1,4} ] * PI_u,'Generator form of [ B(p) ; D(p) ] is OK');

% Minimal generator for PI_u
[S_u,PI_u,iS_u,~] = P_mingen_for_LFR(PI_u.set_vars(p),'lims',p_lims_comp);
[F{1,2},F{1,4},F{4,2},F{4,4}] = pcz_split_matrix([F{1,2},F{1,4};F{4,2},F{4,4}] * S_u, [nx+nu Inf], [nu Inf]);
pcz_lfrzero_report(BD_lfr - [ F{1,2} F{1,4} ] * PI_u,'Minimal generator form of [ B(p) ; D(p) ] is OK');

% Reload PI_2 corresponding to the minimal generator PI_u.
PI_2 = plfr(iS_u(nu+1:end,nu+1:end) * PI_2); % plfr object

% -----------------------

m1 = size(F{1,3},2);
m2 = size(F{1,4},2);
mx = nx + m1;
mu = nu + m2;

% Matrices Fij are already transformed corresponding to the minimal
% generators PI_x and PI_u.
[F{1,:},F{2,:}] = pcz_split_matrix([F{1,:}],[nx ny],[nx nu m1 m2]);


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
pcz_fhzero_report(@(pdp) PI_a_old(pdp)-S_a*PI_a(pdp), pdp, 'Minimal generator for PI_a(p) is OK');

Aa = [
    F{1,1}   F{1,3}   O(nx,m1) O(nx,m1)
    O(m1,nx) O(m1,m1) I(m1)    I(m1)   
    ];

Ca = [ 
    F{2,1}   F{2,3}   O(ny,m1) O(ny,m1)
    ];

Ba = [
    F{1,2}   F{1,4}   O(nx,m1)
    O(m1,nu) O(m1,m2) I(m1)
    ];

Da = [ 
    F{2,2}   F{2,4}   O(ny,m1) 
    ];

Ea = [ I(nx+m1) O(nx+m1,m1+m1) ];

Ga = [ I(nu) O(nu,m2+m1) ];
    
%% Construct annihilators

[N,~,~,N_err] = P_affine_annihilator_for_LFR(PI_x,p_lfr,'lims',p_lims_comp);
pcz_fhzero_report(@(p) N(p)*PI_x(p), p, N_err, 'Annihilator N');

[Na,~,~,Na_err] = P_affine_annihilator_for_LFR(PI_a,[p_lfr;dp_lfr],'lims',pdp_lims_comp);
pcz_fhzero_report(@(pdp) Na(pdp)*PI_a(pdp), pdp, Na_err, 'Annihilator Na');

%% Optimization

% Polytopes (hyperrectangles)
P_v = P_ndnorms_of_X(p_lim);
R_v = P_ndnorms_of_X(dp_lim);
PR_v = P_cartprod(P_v, R_v);

Theta_P = cellfun(@(i){ sdpvar(mx,mx,'symmetric') }, num2cell(0:np));
Q = PGenAffineMatrix(Theta_P,[1;p],p);
dQ = PGenAffineMatrix(Q.get_matrices,[0;dp],pdp);

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
        S_a'*Qa_num*S_a + La*Na_num + Na_num'*La' + eye(size(S_a,2))*1e-10 <= 0 ]; %#ok<AGROW>
    
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
pcz_dispFunction(2,'Model: %s', modelname);
pcz_dispFunction(2,'Overall time: %g', Overall_Time);
pcz_dispFunction(2,'<strong>Gamma = %g</strong> ', gamma);
pcz_dispFunction(2,'size(N)  =%4d, size(PI)   =%4d', size(N.lfrtbx_obj));
pcz_dispFunction(2,'size(Na) =%4d, size(PI_a) =%4d', size(Na.lfrtbx_obj));

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

store_results('Results_All.csv',modelname,0,gamma,sol.solvertime,Overall_Time,sol.info,'Polytopic a. with Finsler')

pcz_dispFunctionEnd(TMP_UFTXCLDbxHBtWRStETWI);

%%
% End of script
pcz_dispFunctionEnd(TMP_BIvQLYBqfFOwenhowBxT);

end