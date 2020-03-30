function method2_descriptor_dual(modelname, A_fh, B_fh, C_fh, D_fh, p_lim, dp_lim, CONDENSED_LMI)
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
% [MAS08] Masubuchi and Suzuki (2008). Gain-scheduled controller synthesis
% based on new LMIs for dissipativity of descriptor LPV systems. IFAC
% Proceedings Volumes, 41(2):9993-9998, 2008.
%
% [MAS99] Masubuchi (1999). An exact solution to parameter-dependent convex
% differential inequalities. 1999 European Control Conference (ECC),
% 1043-1048, 1999.
%
% 1: CONDENSED LMI FORM as proposed by [MAS08]
% 0: Equivalent EXPLICIT LMI FORM as proposed by Author based on  [MAS08]

LMI_form = {'explicit LMI form', 'condensed LMI form'};


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

A = @(x) [ % size: (nx+m1)x(nx+m1)
    F{1,1}              F{1,3}
    diag(dfh(x))*F{3,1} diag(dfh(x))*F{3,3}-I(m1)
    ];

B = @(x) [ % size: (nx+m1)x(nu)
    F{1,2}
    diag(dfh(x))*F{3,2}
    ];

C = [ F{2,1} F{2,3} ]; % size: (nu)x(nx+m1)

D = F{2,2}; % size: (ny)x(nu)

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
% Constraint: [both condensed and explicit]
%
%  Y' E = E' Y >= 0, where E = [ I 0
%                                0 0 ]
%
% Structure for Y:
%
%  T = [ Y1 0  ]
%      [ Y2 Y3 ], where Y1 = Y1'
%
% Y = [ Y11 Y12 ; 0 Y22 ], ahol Y11 = Y11' > 0
Y_Theta = cellfun(@(Y11,Yi2) { [ [ Y11 ; O(m1,nx) ] , Yi2 ] }, ...
    sdpvar(ones(1,nxi)*nx),...
    sdpvar(ones(1,nxi)*(nx+m1),ones(1,nxi)*m1));


%%%
% Constraint: [explicit form]
%
%  E Z1' = 0 and E Z2' = 0
%
% Structure for Z1 and Z2:
%
%  Z1 = [ 0 Z12 ]  and  Z2 = [ 0 Z22 ]
%
% Z1 and Z2 used for the explicit LMI form derived by the author
Z1_Theta = cellfun(@(Z2) { [ O(nu,nx) Z2 ] }, ...
    sdpvar(ones(1,nxi)*nu,ones(1,nxi)*m1,'full') );
Z2_Theta = cellfun(@(Z2) { [ O(ny,nx) Z2 ] }, ...
    sdpvar(ones(1,nxi)*ny,ones(1,nxi)*m1,'full') );


Y = PGenAffineMatrix(Y_Theta,xi,p);
Z1 = PGenAffineMatrix(Z1_Theta,xi,p);
Z2 = PGenAffineMatrix(Z2_Theta,xi,p);

dY = Y.set_channels(dxi,pdp);

gammaSqr = sdpvar;

%%
% (LMI-1) Y(p) E' > 0

P_v = P_ndnorms_of_X(p_lim);
nP = size(P_v,1);

posdef = @(P) P - eye(size(P,1))*1e-8 >= 0;
CONS = [];
for k = 1:nP
    CONS = [ CONS , Y(P_v(k,:)')*E' >= 0 ];
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
Y_eval = Evaluate(Y);
Z1_eval = Evaluate(Z1);
Z2_eval = Evaluate(Z2);

% Get corner index -- mapping from a corner to an index:
% [0 0 0] -> 1
% [0 0 1] -> 2
% [1 1 1] -> 8
get_corner_index = @(corner) corner * (2.^(np-1:-1:0))' + 1;
A_ = @(corner) A_eval( get_corner_index(corner) );
B_ = @(corner) B_eval( get_corner_index(corner) );
Y_ = @(corner) Y_eval( get_corner_index(corner) );
Z1_ = @(corner) Z1_eval( get_corner_index(corner) );
Z2_ = @(corner) Z2_eval( get_corner_index(corner) );


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
    dY_subs_dp = dY.subs(dp,R_v(l,:)',p);

    % Kiertekelem dP(p)-t is, ahogy A(p)-vel es P(p)-vel is tettem.
    dY_eval = Evaluate(dY_subs_dp);

    % Ezekutan dP_-be mar csak a sarokpont indexet kell beirni, hogy
    % megkapjam a mar kiszamolt erteket a megadott indexben.
    dY_ = @(corner) dY_eval( get_corner_index(corner) );

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

        if CONDENSED_LMI

            % Condensed LMI form [MAS08]:
            %
            %  He{ M [ Y' Z'] } + [ -E dY'   0   ]
            %    {   [ 0  I ] }   [   0     Pi2a ]
            %
            TERMS_k = cellfun( @(A,B,Y,dY,Z1,Z2) {
                He{ ...
                    ... Matrix M:
                    [
                        A          , B           , zeros(nx+m1,ny)
                        Pi12'*C    , Pi12'*D     , zeros(nu,ny)
                        Gamma11'*C , Gamma11'*D  , zeros(ny,ny)
                    ] ...
                    * ... Matrix [ Y' Z' ; 0 I ]:
                    [
                        Y'          , Z1'      , Z2'
                        O(nu,nx+m1) , I(nu)    , O(nu,ny)
                        O(ny,nx+m1) , O(ny,nu) , I(ny)
                    ] ...
                } + blkdiag( -E*dY' , Pi22 , -I(ny) )
                ...
            },  A_(fv), B_(fv), Y_(gv), dY_(gv), Z1_(gv), Z2_(gv));

        else

            % Explicit LMI (equivalent form of condensed LMI)
            TERMS_k = cellfun( @(A,B,Y,dY,Z1,Z2) {
                [
                    He(A*Y')-E*dY'      , B+Y*C'*Pi12+A*Z1'             , Y*C'*Gamma11+A*Z2'
                    B'+Pi12'*C*Y'+Z1*A' , Pi22+He(Pi12'*D+Pi12'*C*Z1')  , (D'+Z1*C')*Gamma11+Pi12'*C*Z2'
                    Gamma11'*C*Y'+Z2*A' , Gamma11'*(D+C*Z1')+Z2*C'*Pi12 , He(Gamma11'*C*Z2')-eye(ny)
                ]
                }, A_(fv), B_(fv), Y_(gv), dY_(gv), Z1_(gv), Z2_(gv));

        end

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

Y = value(Y);
Y.name = 'Y';

pcz_2basews(Y);

solver_time = sol.solvertime;
pcz_dispFunction_scalar(gamma, solver_time);
pcz_dispFunction(sol.info);

Y_cell = Y.get_matrices;
Y_names = cellfun(@(i) {sprintf('%s%d',Y.name,i)}, num2cell(0:numel(Y_cell)-1)');
Y_channels = cellfun(@char,num2cell(Y.channels),'UniformOutput',0);

pcz_dispFunction2('%s = %s', Y.name_full, ...
    strjoin(cellfun(@(b,Y) {sprintf('%s * %s',b,Y)}, Y_channels, Y_names), ' + '))

for i = 1:numel(Y_cell)
    pcz_dispFunction_num2str(Y_cell{i}, 'format', '%7.5g','name',Y_names{i})
end


store_results('Descriptor_Results.csv', modelname, 0, gamma, sol.solvertime, Overall_Time, ...
    sol.info, [ 'dual ' LMI_form{CONDENSED_LMI+1} ])

store_results('Results_All.csv', modelname, 0, gamma, sol.solvertime, Overall_Time, ...
    sol.info, [ 'Descriptor - dual ' LMI_form{CONDENSED_LMI+1} ])

%%

pcz_dispFunctionEnd(TMP_yMJTNGorRueTwiaHttwA);

end

