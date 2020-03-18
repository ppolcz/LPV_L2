function LPV_L2anal_Masubuchi_primal(A_fh, B_fh, C, D, K, Use_MinLFR, x_lim, p_expr, p_lim, dp_lim, varargin)
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
% [MAS99] Masubuchi (1999). An exact solution to parameter-dependent convex
% differential inequalities. 1999 European Control Conference (ECC),
% 1043-1048, 1999. 
% 
% This implementation is an extension of [MAS03, Lemma 1] for linear
% parameter varying systems.


%%

I = @(varargin) eye(varargin{:});
O = @(varargin) zeros(varargin{:});

He = he;

TMP_yMJTNGorRueTwiaHttwA = pcz_dispFunctionName;

%% Build up a descriptor model representation

% Find out np, nw and nz
[nz,nx] = size(C);
[~,nw] = size(D);
np = size(p_lim,1);

P_generate_symvars(nx,np,nw,nz);

[p_lfr,p_lfr_cell] = pcz_generateLFRStateVector('p',p_lim);

M_fh = @(varargin) [ 
    A_fh(varargin{:}) - B_fh(varargin{:})*K , B_fh(varargin{:})
    C                                       , D
    ];
M_lfr = M_fh(p_lfr_cell{:});

if Use_MinLFR == 1
    
    % Here, we use the numerical n-D order reduction proposed by
    % 
    % [DAK97] D'Andrea and Khatri (1997). Kalman decomposition of linear
    % fractional transformation representations and minimality. Proceedings
    % of the 1997 American Control Conference (Cat. No.97CH36041),
    % 6:3557-3561 vol.6, 1997.

    M_lfr = minlfr(M_lfr);

elseif Use_MinLFR == 2
    
    % Due to the nice properties of the model, we can also use the
    % symbolical LFR transformation of
    % 
    % [PPS18a] Polcz, Péni and Szederkényi (2018). Computational method for
    % estimating the domain of attraction of discrete-time uncertain
    % rational systems. European Journal of Control, , 2018.
    %     
    % [PPS18b] Polcz, Péni and Szederkényi (2018). Reduced linear
    % fractional representation of nonlinear systems for stability
    % analysis. IFAC-PapersOnLine, 51(2):37-42, 2018.

    M_plfr_red = P_LFR_reduction(M_lfr,[x;w]);
    l = M_plfr_red.LFR_reduced;
    M_lfr = plfr(l.A,l.B,l.C,l.D,l.Delta,p_lim);
    clear l
end

M_plfr = plfr(M_lfr);

pcz_lfrzero_report(M_plfr - M_fh(p_lfr_cell{:}), ...
    'The reduced LFR corresponds to the original one')

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
[F{:}] = pcz_split_matrix(M_plfr.M,[nx,nz,m1],[nx,nw,m1],'RowWise',0);
Delta = M_plfr.Delta;

%%% 
% Descriptor model:
% 
%   E dX = A(p) Xi + B(p) u,     where X = [ x
%      y = C    Xi + D    u,                 pi ] in R^n
% 
E = blkdiag(I(nx), O(m1));

A_sym = [
    F{1,1}       F{1,3}
    Delta*F{3,1} Delta*F{3,3}-I(m1)
    ];
n = nx + m1;

B_sym = [
    F{1,2}
    Delta*F{3,2}
    ];

C = [ F{2,1} F{2,3} ];

D = F{2,2};

%%

delta = diag(Delta);
dfh = matlabFunction(delta,'vars',{p});

A = @(x) [
    F{1,1}              F{1,3}
    diag(dfh(x))*F{3,1} diag(dfh(x))*F{3,3}-I(m1)
    ];
    
B = @(x) [
    F{1,2}
    diag(dfh(x))*F{3,2}
    ];


%%

% Generate every multi-affine monoms of p
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
%  P' E = E' P >= 0, where E = [ I 0
%                                0 0 ]
% 
% Structure for P:
% 
%  P = [ P1 0
%        P2 P3 ], where P1 = P1'
%     
P_Theta = cellfun(@(P1,P23) { [ P1 O(nx,m1) ; P23 ] }, ...
    sdpvar(ones(1,nxi)*nx,ones(1,nxi)*nx,'symmetric'),...
    sdpvar(ones(1,nxi)*m1,ones(1,nxi)*(nx+m1),'full'));

%%%
% Constraint:
% 
%  Z1' E = 0 and Z2' E = 0
% 
% Structure for Z1 and Z2:
% 
%  Z1 = [ 0       and  Z2 = [ 0
%         Z12 ],              Z22 ]
%     
Z1_Theta = cellfun(@(Z) { [ O(nx,nw) ; Z ] }, ...
    sdpvar(ones(1,nxi)*m1,ones(1,nxi)*nw,'full') );
Z2_Theta = cellfun(@(Z) { [ O(nx,nz) ; Z ] }, ...
    sdpvar(ones(1,nxi)*m1,ones(1,nxi)*nz,'full') );

P = PGenAffineMatrix(P_Theta,xi,p);
Z1 = PGenAffineMatrix(Z1_Theta,xi,p);
Z2 = PGenAffineMatrix(Z2_Theta,xi,p);

dP = P.set_channels(dxi,pdp);

gammaSqr = sdpvar;

%%
%(LMI-1) P(p) > 0

P_v = P_ndnorms_of_X(p_lim);
nP = size(P_v,1);

posdef = @(P) P - eye(size(P,1))*1e-3 >= 0;
CONS = [];
for k = 1:nP
    CONS = [ CONS , E'*P(P_v(k,:)') >= 0 ];
end

%%
% (LMI-2) He{ P(p) A(p) } + dP(p,dp) < 0
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
P_eval = Evaluate(P);
Z1_eval = Evaluate(Z1);
Z2_eval = Evaluate(Z2);

% Get corner index -- mapping from a corner to an index:
% [0 0 0] -> 1
% [0 0 1] -> 2
% [1 1 1] -> 8
get_corner_index = @(corner) corner * (2.^(np-1:-1:0))' + 1;
A_ = @(corner) A_eval( get_corner_index(corner) );
B_ = @(corner) B_eval( get_corner_index(corner) );
P_ = @(corner) P_eval( get_corner_index(corner) );
Z1_ = @(corner) Z1_eval( get_corner_index(corner) );
Z2_ = @(corner) Z2_eval( get_corner_index(corner) );


% I render the corner combinations of [MAS99] in advance:
combs = repmat({ { {0 0} , {1 1} , num2cell([ 0 1 ; 1 0 ],1) } },[1 np]);
combs = allcomb(combs{:});


% Matrices of the supply function 
Pi22 = -eye(nw)*gammaSqr;
Pi12 = 0;
Gamma11 = eye(nz);

stringify1 = @(fv) sprintf('(%s)', strjoin(cellfun(@(n) {num2str(n)}, num2cell(fv)),','));
stringify = @(fv) cellfun(@(c) {stringify1(c)}, num2cell(fv,2));

% I walk through the corner points of R
for l = 1:Nr_Combs_dp

    pcz_dispFunction2('%d. corner %s of R\n', l, stringify1(R_v_corners(l,:)))

    % Behelyettesitem dP(p,dp=vR) = dP(p), ahol vR az R egy sarokpontja.
    % Ezutan dP(p) mar csak p-tol fog fuggeni.
    dP_subs_dp = dP.subs(dp,R_v(l,:)',p);
    
    % Kiertekelem dP(p)-t is, ahogy A(p)-vel es P(p)-vel is tettem.
    dP_eval = Evaluate(dP_subs_dp);
    
    % Ezekutan dP_-be mar csak a sarokpont indexet kell beirni, hogy
    % megkapjam a mar kiszamolt erteket a megadott indexben.
    dP_ = @(corner) dP_eval( get_corner_index(corner) );

    for k = 1:Nr_Combs_p

        combk = vertcat(combs{k,:});

        fv = allcomb(combk{:,1});
        gv = allcomb(combk{:,2});

        pcz_dispFunction2('    %d.%d. Combination nr. %d:', l,k,k)

        str = cellfun(@(f,g) { [f ' and ' g] }, stringify(fv),stringify(gv));
        for j = 1:numel(str)
            pcz_dispFunction2('        %d.%d.%d. Corner points of term %2d: %s', l,k,j, j, str{j})
        end
        pcz_dispFunction2 ' '

%         TERMS_k = cellfun( @(A,B,P,dP) { 
%             [ 
%             He( P'*A ) + dP'*E , P'*B + C'*Pi12     , C'*Gamma11
%             B'*P + Pi12'*C     , Pi22 + He(D'*Pi12) , D'*Gamma11
%             Gamma11'*C         , Gamma11'*D         , -eye(nz)
%             ]
%             }, A_(fv), B_(fv), P_(gv), dP_(gv));

        TERMS_k = cellfun( @(A,B,P,dP,Z1,Z2) { 
            [ 
            He(P'*A)+dP'*E     , P'*B+C'*Pi12+A'*Z1     , C'*Gamma11+A'*Z2
            B'*P+Pi12'*C+Z1'*A , Pi22+He(D'*Pi12+Z1'*B) , D'*Gamma11+B'*Z2
            Gamma11'*C+Z2'*A   , Gamma11'*D+Z2'*B       , -eye(nz)
            ]
            }, A_(fv), B_(fv), P_(gv), dP_(gv), Z1_(gv), Z2_(gv));

        if numel(TERMS_k) == 1
            CONS = [ CONS , posdef(-TERMS_k{1}) ];
        else
            CONS = [ CONS , posdef(-sum( cat(3,TERMS_k{:}), 3 )) ];
        end

    end

end

CONS

sol = optimize(CONS,gammaSqr)
gamma = sqrt(value(gammaSqr))

P = value(P)

store_results('Descriptor_Results.csv', x_lim, gamma, 0, sol.solvertime, ...
    sol.info, 'primal', Use_MinLFR)

return

dP_subs_dp = value(dP_subs_dp);

eigP = @(p) min(eig(P(p))); 

eigAP = @(p) max(eig(He{ P(p)*A(p) } + dP_subs_dp(p)));

% ELLENORZES

resolution = {
    10
    10
    10
    };
p_lim_cell = num2cell(num2cell(p_lim),2);

p_lspace = cellfun(@(res,lims) {linspace(lims{:},res)}, resolution, p_lim_cell);

p_mesh = cell(size(p_lspace));

[p_mesh{:}] = ndgrid(p_lspace{:});

p_mesh = cellfun(@(a) {a(:)'}, p_mesh);
pp = num2cell(vertcat(p_mesh{:}),1);
eigP_num = cellfun(eigP, pp);
eigAP_num = cellfun(eigAP, pp);

figure, hold on
plot(eigP_num)
plot(eigAP_num)

legend 'min eig of P' 'max eig of He\{AP\}+dP'


%%

pcz_dispFunctionEnd(TMP_yMJTNGorRueTwiaHttwA);

end

function cross_corner_evaluation
%%

    
end
