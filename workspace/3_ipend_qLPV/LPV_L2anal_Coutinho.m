%% LPV_L2anal_Coutinho
%
%  File: LPV_L2anal_Coutinho.m
%  Directory: 1_PhD_projects/22_Hinf_norm/LPV_TAC_v3_newmod
%  Author: Peter Polcz (ppolcz@gmail.com)
%
%  Created on 2019. September 16. (2019a)
%

%%
% Automatically generated stuff

global SCOPE_DEPTH VERBOSE LATEX_EQNR
SCOPE_DEPTH = 0;
VERBOSE = 1;
LATEX_EQNR = 0;


try c = evalin('caller','persist'); catch; c = []; end
persist = Persist(mfilename('fullpath'), c); clear c;
persist.backup();
%clear persist

disp 'EBBEN A FORMABAN NEM ALKALMAZHATO LPV RENDSZERRE!'

return

%%

I = @(varargin) eye(varargin{:});
O = @(varargin) zeros(varargin{:});

He = he;

Use_MinLFR = 0;

%%

P_generate_symvars_v10(4,3,2,2);

p_lim = [
    -1 2
    -1 2
    0 2
    ];

dp_lim = [
    -10 10
    -1 1
    -5 5
    ];

[p_lfr,p_lfr_cell] = pcz_generateLFRStateVector('p',p_lim);

M_fh = @(p1,p2,p3) [
    -3+p1       3+p1/(p3^2 + 0.5*p1 + 1)     0.1*p3   p3^2*p2       , 0      0
    0          -1-p2^2                       5        0             , 1+p2^2 0
    -1/(5-p2)  0                             -4+p1    0             , 0      0
    0          0.1                           0        -5+1/(p1 + 2) , 0      2+p1/(p1 + 2)
    ............................................................... , ...................
    1/(5-p2)   0                             0        0             , 0      0
    0          0                             p1+1     0             , 0      0
    ];

M_sym = M_fh(p_cell{:});

M_lfr = M_fh(p_lfr_cell{:});

if Use_MinLFR
    M_lfr = minlfr(M_lfr);
end

M_plfr = plfr(M_lfr);
m1 = M_plfr.m1;

% El Ghaoui es Scorletti jeloleseivel
[A,Bu,Bp , Cy,Dyu,Dyp , Cq,Dqu,Dqp] = pcz_split_matrix(M_plfr.M,[nx,nz,m1],[nx,nw,m1]);
Delta = M_plfr.Delta;

delta_fh = matlabFunction(diag(Delta),'vars',{p});

%% TOROLHETO

% F = cell(3,3);
% [F{:}] = pcz_split_matrix(M_plfr.M,[nx,nz,m1],[nx,nw,m1],'RowWise',0);
% m1 = size(Delta,1);

% [A1,A3,A2] = deal(F{1,1:3});
% [E1,E3,E2] = deal(F{2,1:3});


% PI_1 = @(x) diag(delta_fh(x))*F{3,1};
% PI_1 = @(x) diag(delta_fh(x))*F{3,3} - I(m1);
% PI_1 = @(x) diag(delta_fh(x))*F{3,2};


%% Model szimbolikus ellenorzese

%{

ABCD_sym = [ A Bu ; Cy Dyu ] + [ Bp ; Dyp ] / (I(m1) - Delta*Dqp) * Delta * [ Cq Dqu ];

pcz_symzero_report(simplify(ABCD_sym - M_sym), 'Az LFR modell rendben van')


%}

%% Coutinho-fele DAR

A1 = A;
A2 = Bp;
A3 = Bu;

E1 = Cy;
E2 = Dyp;
E3 = Dyu;

PI_1 = @(p) diag(delta_fh(p))*Cq;
PI_2 = @(p) diag(delta_fh(p))*Dqp - I(m1);
PI_3 = @(p) diag(delta_fh(p))*Dqu;

PI_1_sym = Delta*Cq;
PI_2_sym = Delta*Dqp - I(m1);
PI_3_sym = Delta*Dqu;

Pi = - PI_2_sym \ (PI_1_sym * x + PI_3_sym * w);

%% Storage function felépítése

Theta = kron(p,I(nx));
n_Theta = size(Theta,1);

zeta = [
    Theta * x
    x
    ];


% zeta affin annihilatora
Psi_1 = [ I(n_Theta) -Theta ];

ZEROS = simplify(Psi_1 * zeta)

% N1*zeta == x
N1 = [ O(nx,n_Theta) I(nx) ];

ZEROS = N1*zeta - x

% x affin annihilatora
Nx_fh = P_affine_annihilator(x,x,x,'sym',1);
Nx = sym(Nx_fh);
n_Nx = size(Nx,1);

ZERO = Nx * x


% col[1,zeta] affin annihilatora
% Ez csak a boundary condition-hoz kell (ami most nincs)
% Psi_2 = [
%     x            -N1
%     O(n_Theta,1) Psi_1
%     ];
% 
% ZERO = Psi_2 * [1 ; zeta]

% col[zeta,Pi,w] affin annihilatora
Psi_3_fh = @(p) [
    PI_1(p)*N1       PI_2(p)        PI_3(p)
    ];

% ZERO = simplify(Psi_3 * [zeta;Pi;w])

%% Storage function deriváltja

Theta_tilde = [
    Theta + kron(I(nx),p)
    I(nx)
    ];

%{

P = pcz_sym_symmetric('p',nx+n_Theta);


V_elm = zeta'*P*zeta;

dV_elm = jacobian(V_elm,x)*dx;

dV_eq = 2*zeta'*P*Theta_tilde*dx;

ZERO = simplify(dV_elm - dV_eq);

%}

%% Polytope

% outer X polytope: az alaphalmaz.
% v   : vertices
% fci : facet corner indices
% ak  : vectors from the origin, perpendicular to facets
[P_v, P_ak, P_fci] = P_ndnorms_of_X(x_lim);
P_NrF = size(P_fci,1);
P_NrC = size(P_fci,2);

%%% Optimization

sdpvar gamma

P = sdpvar(nx+n_Theta);

L = sdpvar(size(Psi_1,2),size(Psi_1,1));

Psi_3_dummy = Psi_3_fh(double(x*0));
W = sdpvar(size(Psi_3_dummy,2),size(Psi_3_dummy,1));

CONS = [ mu >= 0 ];

Psi_1_fh = matlabFunction(Psi_1,'vars',{p});
% Psi_3_fh = matlabFunction(Psi_3,'vars',{x});

Theta_tilde_fh = matlabFunction(Theta_tilde,'vars',{x});

for i = 1:size(P_v,1)
    Psi_1_num = Psi_1_fh(P_v(i,:)');
    Psi_3_num = Psi_3_fh(P_v(i,:)');

    Theta_tilde_num = Theta_tilde_fh(P_v(i,:)');

    CONS = [CONS ,

        P + L*Psi_1_num + Psi_1_num'*L' - eye(size(P))*1e-10 >= 0

        [
            He(P*Theta_tilde_num*A1*N1) P*Theta_tilde_num*A2 P*Theta_tilde_num*A3 N1'*E1'
            A2'*Theta_tilde_num'*P      O(m1)                O(m1,nw)             E2'
            A3'*Theta_tilde_num'*P      O(nw,m1)             -gamma*I(nw)         E3'
            E1*N1                       E2                   E3                   -gamma*I(nz)
        ] + He( [ W ; O(size(Psi_3_dummy,1),nz)' ] * [ Psi_3_num O(size(Psi_3_num,1),nz) ] ) <= 0

        ];

end

sdps = sdpsettings('solver', 'mosek', 'verbose', true, 'cachesolvers', true);
sol = optimize(CONS, gamma, sdps)

double(gamma)

%%
persist.stoplog;
