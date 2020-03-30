function method2_descriptor_primal(modelname, A_fh, B_fh, C_fh, D_fh, p_lim, dp_lim, varargin)
%%
%  
%  File: LPV_L2anal_Masubuchi_dual.m
%  Directory: 1_PhD_projects/22_Hinf_norm/qLPV_ipend
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2020. March 13. (2019b)
%  Minor review on 2020. March 17. (2019b)
%
% 
% Based on:
% 
% [MAS03] Masubuchi, Akiyama and Saeki (2003). Synthesis of output feedback
% gain-scheduling controllers based on descriptor LPV system
% representation. 42th IEEE Conference on Decision and Control, 6:6115-6120
% Vol.6, 2003
% 
% [MAS06] Masubuchi (2006). Dissipativity inequalities for continuous-time
% descriptor systems with applications to synthesis of control gains.
% Systems & Control Letters, 55(2):158-164, 2006.
% 
% [MAS08] Masubuchi and Suzuki (2008). Gain-scheduled controller synthesis
% based on new LMIs for dissipativity of descriptor LPV systems. IFAC
% Proceedings Volumes, 41(2):9993-9998, 2008.
% 
% [MAS99] Masubuchi (1999). An exact solution to parameter-dependent convex
% differential inequalities. 1999 European Control Conference (ECC),
% 1043-1048, 1999. 
% 
% This implementation is the equivalent ``primal'' version of the condensed
% dual LMI formulation of [MAS08]


%%

I = @(varargin) eye(varargin{:});
O = @(varargin) zeros(varargin{:});

He = he;

TMP_yMJTNGorRueTwiaHttwA = pcz_dispFunctionName;

%% Build up a descriptor model representation

[~,B_lfr,C_lfr,~,~,~,M_lfr] = helper_convert(A_fh, B_fh, C_fh, D_fh, p_lim);

% Dimensions
nx = size(B_lfr,1);
nu = size(B_lfr,2);
ny = size(C_lfr,1);
np = size(p_lim,1);

% Generate symbolic variables
P_generate_symvars(nx,np,nu,ny);

% Here, we use the numerical n-D order reduction proposed by
%
% [DAK97] D'Andrea and Khatri (1997). Kalman decomposition of linear
% fractional transformation representations and minimality. Proceedings of
% the 1997 American Control Conference (Cat. No.97CH36041), 6:3557-3561
% vol.6, 1997.
M_plfr = plfr(minlfr(M_lfr));
m1 = M_plfr.m1;

%%%
% LFR model:
%
%   dx = F11 x + F12 u + F13 pi
%    y = F21 x + F22 u + F23 pi
%  eta = F31 x + F32 u + F33 pi
%   pi = Delta eta
%
F = cell(3,3);
[F{:}] = pcz_split_matrix(M_plfr.M,[nx,ny,m1],[nx,nu,m1],'RowWise',0);
Delta = M_plfr.Delta;

% Helper
delta = diag(Delta);
dfh = matlabFunction(delta,'vars',{p});

%%
% Descriptor model:
%
%   E dX = A(p) Xi + B(p) u,     where X = [ x
%      y = C    Xi + D    u,                 pi ] in R^n
%
E = blkdiag(I(nx), O(m1));

A = @(x) [
    F{1,1}              F{1,3}
    diag(dfh(x))*F{3,1} diag(dfh(x))*F{3,3}-I(m1)
    ];

B = @(x) [
    F{1,2}
    diag(dfh(x))*F{3,2}
    ];

C = [ F{2,1} F{2,3} ];

D = F{2,2};

%% Render multiaffine monoms and multiaffine decision variables
%
%  degree 1 monoms: p1, p2, ...
%  degree 2 monoms: p1*p2, p1*p3, p1*p4, p2*p3, ....
%  degree 3 monoms: p1*p2*p3, p1*p2*p4, ....
monoms = cellfun( @(k) { prod(combnk(p,k),2) }, num2cell(1:np) );

xi = [
    1
    vertcat(monoms{:})
    ];

dxi = jacobian(xi,p)*dp;

nxi = 2^np;

%%%
% Constraint:
% 
%  X' E = E' X >= 0, where E = [ I 0
%                                0 0 ]
% 
% Structure for X:
% 
%  X = [ X1 0
%        X2 X3 ], where X1 = X1'
%
X_Theta = cellfun(@(X1,X23) { [ X1 O(nx,m1) ; X23 ] }, ...
    sdpvar(ones(1,nxi)*nx),...
    sdpvar(ones(1,nxi)*m1,ones(1,nxi)*(nx+m1),'full'));

%%%
% Constraint:
%
%  W1' E = 0 and W2' E = 0
%
% Structure for W1 and W2:
%
%  W1 = [ 0       and  W2 = [ 0
%         W12 ],              W22 ]
%
W1_Theta = cellfun(@(W) { [ O(nx,nu) ; W ] }, ...
    sdpvar(ones(1,nxi)*m1,ones(1,nxi)*nu,'full') );
W2_Theta = cellfun(@(W) { [ O(nx,ny) ; W ] }, ...
    sdpvar(ones(1,nxi)*m1,ones(1,nxi)*ny,'full') );

X = PGenAffineMatrix(X_Theta,xi,p);
W1 = PGenAffineMatrix(W1_Theta,xi,p);
W2 = PGenAffineMatrix(W2_Theta,xi,p);

dX = X.set_channels(dxi,pdp);

gammaSqr = sdpvar;

%%
% (LMI-1) E' X(p) > 0

P_v = P_ndnorms_of_X(p_lim);
nP = size(P_v,1);

posdef = @(P) P - eye(size(P,1))*1e-8 >= 0;
CONS = [];
for k = 1:nP
    CONS = [ CONS , E'*X(P_v(k,:)') >= 0 ];
end

%%
% (LMI-2) He{ X(p) A(p) } + E' dX(p,dp) < 0
%
% Implementation of the cross-corner method of [MAS99].
%
% #Cross-corner method with respect to p.
% #Simply corner points with respect to dp.

R_v = P_ndnorms_of_X(dp_lim);
R_v_corners = P_ndnorms_of_X(ones(size(dp_lim)) * [0 0;0 1]);
nR = size(R_v,1);

% Number of combinations
Nr_Combs_p = 3^np;
Nr_Combs_dp = 2^np;

% Evaluation in each corner point
Evaluate = @(f) cellfun( @(x_num) { f(x_num) }, num2cell(P_v',1));

A_eval = Evaluate(A);
B_eval = Evaluate(B);
X_eval = Evaluate(X);
W1_eval = Evaluate(W1);
W2_eval = Evaluate(W2);

% Get corner index -- mapping from a corner to an index:
% [0 0 0] -> 1
% [0 0 1] -> 2
% [1 1 1] -> 8
get_corner_index = @(corner) corner * (2.^(np-1:-1:0))' + 1;
A_ = @(corner) A_eval( get_corner_index(corner) );
B_ = @(corner) B_eval( get_corner_index(corner) );
X_ = @(corner) X_eval( get_corner_index(corner) );
W1_ = @(corner) W1_eval( get_corner_index(corner) );
W2_ = @(corner) W2_eval( get_corner_index(corner) );


% I render the corner combinations of [MAS99] in advance:
combs = repmat({ { {0 0} , {1 1} , num2cell([ 0 1 ; 1 0 ],1) } },[1 np]);
combs = allcomb(combs{:});


% Matrices related to the supply function
Pi22 = -eye(nu)*gammaSqr;
Pi12 = O(ny,nu);
Gamma11 = eye(ny);

stringify1 = @(fv) sprintf('(%s)', strjoin(cellfun(@(n) {num2str(n)}, num2cell(fv)),','));
stringify = @(fv) cellfun(@(c) {stringify1(c)}, num2cell(fv,2));

% I walk through the corner points of R
for l = 1:Nr_Combs_dp

    % pcz_dispFunction2('%d. corner %s of R\n', l, stringify1(R_v_corners(l,:)))

    % Behelyettesitem dP(p,dp=vR) = dP(p), ahol vR az R egy sarokpontja.
    % Ezutan dP(p) mar csak p-tol fog fuggeni.
    dX_subs_dp = dX.subs(dp,R_v(l,:)',p);

    % Kiertekelem dP(p)-t is, ahogy A(p)-vel es P(p)-vel is tettem.
    dX_eval = Evaluate(dX_subs_dp);

    % Ezekutan dP_-be mar csak a sarokpont indexet kell beirni, hogy
    % megkapjam a mar kiszamolt erteket a megadott indexben.
    dX_ = @(corner) dX_eval( get_corner_index(corner) );

    for k = 1:Nr_Combs_p

        combk = vertcat(combs{k,:});

        fv = allcomb(combk{:,1});
        gv = allcomb(combk{:,2});

        % pcz_dispFunction2('    %d.%d. Combination nr. %d:', l,k,k)

        str = cellfun(@(f,g) { [f ' and ' g] }, stringify(fv),stringify(gv));
        
        % for j = 1:numel(str)
        %     pcz_dispFunction2('        %d.%d.%d. Corner points of term %2d: %s', l,k,j, j, str{j})
        % end
        % pcz_dispFunction2 ' '

        TERMS_k = cellfun( @(A,B,X,dX,W1,W2) { 
            [
            He(X'*A)+dX'*E     , X'*B+C'*Pi12+A'*W1     , C'*Gamma11+A'*W2
            B'*X+Pi12'*C+W1'*A , Pi22+He(D'*Pi12+W1'*B) , D'*Gamma11+B'*W2
            Gamma11'*C+W2'*A   , Gamma11'*D+W2'*B       , -eye(ny)
            ]
            }, A_(fv), B_(fv), X_(gv), dX_(gv), W1_(gv), W2_(gv));

        if numel(TERMS_k) == 1
            CONS = [ CONS , posdef(-TERMS_k{1}) ];
        else
            CONS = [ CONS , posdef(-sum( cat(3,TERMS_k{:}), 3 )) ];
        end

    end

end

sol = optimize(CONS,gammaSqr,sdpsettings('solver','mosek'));
gamma = sqrt(value(gammaSqr));

Overall_Time = toc(TMP_yMJTNGorRueTwiaHttwA);

X = value(X);
X.name = 'X';

pcz_2basews(X);

solver_time = sol.solvertime;
pcz_dispFunction_scalar(gamma, solver_time);
pcz_dispFunction(sol.info);

X_cell = X.get_matrices;
X_names = cellfun(@(i) {sprintf('%s%d',X.name,i)}, num2cell(0:numel(X_cell)-1)');
X_channels = cellfun(@char,num2cell(X.channels),'UniformOutput',0);

pcz_dispFunction2('%s = %s', X.name_full, ...
    strjoin(cellfun(@(b,X) {sprintf('%s * %s',b,X)}, X_channels, X_names), ' + '))

for i = 1:numel(X_cell)
    pcz_dispFunction_num2str(X_cell{i}, 'format', '%7.5g','name',X_names{i})    
end


store_results('Descriptor_Results.csv', modelname, 0, gamma, sol.solvertime, Overall_Time, ...
    sol.info, 'primal')

store_results('Results_All.csv', modelname, 0, gamma, sol.solvertime, Overall_Time, ...
    sol.info, 'Descriptor - primal')

%%

pcz_dispFunctionEnd(TMP_yMJTNGorRueTwiaHttwA);

end

