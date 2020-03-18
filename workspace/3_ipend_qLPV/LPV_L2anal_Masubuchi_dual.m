function LPV_L2anal_Masubuchi_dual(A_fh, B_fh, C, D, K, Use_MinLFR, x_lim, p_expr, p_lim, dp_lim, varargin)
%% LPV_L2anal_Masubuchi
%  
%  File: LPV_L2anal_Masubuchi_dual.m
%  Directory: 1_PhD_projects/22_Hinf_norm/qLPV_ipend
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2020. March 13. (2019b)
%

%% Eredmenyek
% 
% 2019.09.09. (szeptember  9, hétfő), 04:52
% 
% Nincs minlfr, nincs Z1, Z2: GAMMA = 2.1289, T = 1339 sec
% Van minlfr, nincs Z1, Z2: Gamma = 2.1268, T = 170 sec
% Van minlfr, van Z1, Z2: Gamma = 2.1247, T = 255 sec


%%

I = @(varargin) eye(varargin{:});
O = @(varargin) zeros(varargin{:});

He = he;

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

%% Modell felepitese

% Find out np, nw and nz
[nz,nx] = size(C);
[~,nw] = size(D);
np = size(p_lim,1);

P_generate_symvars_v10(nx,np,nw,nz);

[p_lfr,p_lfr_cell] = pcz_generateLFRStateVector('p',p_lim);

% LFR realization of matrix B(p):
B_lfr = B_fh(p_lfr_cell{:});

% LFR realization of matrix A(p) - B(p)*K:
Ak_lfr = A_fh(p_lfr_cell{:}) - B_lfr*K;

M_lfr = [
    Ak_lfr B_lfr
    C      D
    ];

if Use_MinLFR
    M_lfr = minlfr(M_lfr);
end

M_plfr = plfr(M_lfr);

m1 = M_plfr.m1;

F = cell(3,3);
[F{:}] = pcz_split_matrix(M_plfr.M,[nx,nz,m1],[nx,nw,m1],'RowWise',0);
Delta = M_plfr.Delta;

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


%% Disszipacios betaplalasi fuggveny matrixai
% 
%  Az ami a cikkben m1 = dim(w) = nw
%                   p1 = dim(z) = nz  

gammaSqr = sdpvar;

Pi22 = -eye(nw)*gammaSqr;
Pi12 = O(nz,nw);
Gamma11 = eye(nz);
q = nz;

Pi1s = [ Pi12 Gamma11 ];
Psi = [ I(nw) zeros(nw,q) ];
Pi2a = blkdiag(Pi22,-I(q));

%% [Mashubuchi and Suzuki, 2008, Eq.(4) es (5) kozott (Acl tilde...)]

% Ezt csak azert, hogy szimbolikusan lassuk mink van
M_sym = [
    A_sym   , B_sym*Psi
    Pi1s'*C , Pi1s'*D*Psi
    ];

M = @(p) [
    A(p)    , B(p)*Psi
    Pi1s'*C , Pi1s'*D*Psi
    ];

%% Multiaffin monomok deklaralasa

monoms = cellfun( @(k) { prod(combnk(p,k),2) }, num2cell(1:np) );

xi = [
    1
    vertcat(monoms{:})
    ];

dxi = jacobian(xi,p)*dp;

nxi = 2^np;

%% Szabad dontesi valtozok letrehozasa (Y es Z)

% Y = [ Y11 Y12 ; 0 Y22 ], ahol Y11 = Y11' > 0
Y_Theta = cellfun(@(Y11,Yi2) { [ [ Y11 ; O(m1,nx) ] , Yi2 ] }, ...
    sdpvar(ones(1,nxi)*nx),...
    sdpvar(ones(1,nxi)*(nx+m1),ones(1,nxi)*m1));

% Z = [ 0 Z2 ] 
Z_Theta = cellfun(@(Z) { [ O(nw+q,nx) Z ] }, ...
    sdpvar(ones(1,nxi)*(nw+q),ones(1,nxi)*m1) );

Y = PGenAffineMatrix(Y_Theta,xi,p);
Z = PGenAffineMatrix(Z_Theta,xi,p);

dY = Y.set_channels(dxi,pdp);


%%
% (LMI-1) Y(p) E' > 0

P_v = P_ndnorms_of_X(p_lim);
nP = size(P_v,1);

posdef = @(P) P - eye(size(P,1))*1e-3 >= 0;
CONS = [];
for k = 1:nP
    CONS = [ CONS , Y(P_v(k,:)')*E' >= 0 ];
end

%%
% (LMI-2) He{ P(p) A(p) } + dP(p,dp) < 0
% 
% p-szerint Masubuchi keresztsarokpontos modszere, dp-szerint sarokpontok

R_v = P_ndnorms_of_X(dp_lim);
nR = size(R_v,1);

% Kombinaciok szama
Nr_Combs_p = 3^np;
Nr_Combs_dp = 2^np;

% Kiertekeles minden sarokpontban
Kiertekeles = @(f) cellfun( @(x_num) { f(x_num) }, num2cell(P_v',1));

M_eval = Kiertekeles(M);
Y_eval = Kiertekeles(Y);
Z_eval = Kiertekeles(Z);

get_corner_index = @(corner) corner * (2.^(np-1:-1:0))' + 1;
M_ = @(corner) M_eval( get_corner_index(corner) );
Y_ = @(corner) Y_eval( get_corner_index(corner) );
Z_ = @(corner) Z_eval( get_corner_index(corner) );


% Elore legyartom [Masubuchi, 1999] kombinacioit
combs = repmat({ { {0 0} , {1 1} , num2cell([ 0 1 ; 1 0 ],1) } },[1 np]);
combs = allcomb(combs{:});


% Vegigmegyek R sarokpontjain
for l = 1:Nr_Combs_dp

    % Behelyettesitem dP(p,dp=vR) = dP(p), ahol vR az R egy sarokpontja.
    % Ezutan dP(p) mar csak p-tol fog fuggeni.
    dY_subs_dp = dY.subs(dp,R_v(l,:)',p);
    
    % Kiertekelem dP(p)-t is, ahogy A(p)-vel es P(p)-vel is tettem.
    dY_eval = Kiertekeles(dY_subs_dp);
    
    % Ezekutan dP_-be mar csak a sarokpont indexet kell beirni, hogy
    % megkapjam a mar kiszamolt erteket a megadott indexben.
    dY_ = @(corner) dY_eval( get_corner_index(corner) );

    fprintf('Kiertekeles a `jellegzetes` pontokban:\n\n')
    stringify = @(fv) cellfun(@(s) {[ '(' strrep(s,'  ',',') ')' ]}, num2cell(num2str(fv),2));
    for k = 1:Nr_Combs_p

        combk = vertcat(combs{k,:});

        fv = allcomb(combk{:,1});
        gv = allcomb(combk{:,2});

        fprintf('%d. kombinacio:\n\n', k)

        str = cellfun(@(f,g) { [f ' es ' g] }, stringify(fv),stringify(gv));
        for j = 1:numel(str)
            fprintf('    %2d. tag sarokponjai: %s\n', j, str{j})
        end
        disp ' '

        TERMS_k = cellfun( @(M,Y,dY,Z) { 
            He{ M * [ Y' Z' ; O(nw+q,nx+m1) I(nw+q) ] } + blkdiag( -E*dY' , Pi2a )
            }, M_(fv), Y_(gv), dY_(gv), Z_(gv));

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

Y = value(Y)


%% Store results

s.Stamp = persist.stamp;
s.x1_max = x_lim(1,2);
s.x2_max = x_lim(2,2);
s.x3_max = x_lim(3,2);
s.Gamma = gamma;
s.Solver_Time = sol.solvertime;
s.Solver_info = sol.info;
s.Method = 'dual';
s.Use_MinLFR = Use_MinLFR;

Results_csv = persist.file_simple('txt','Mashubuchi_Results.csv');
pcz_dispFunction('Mashubuchi_Results.csv: `%s''', Results_csv);

if exist(Results_csv, 'file')
    Results = readtable(Results_csv,'Delimiter','|');
end

if ~exist('Results', 'var')
    Cols = { 'Stamp', 'x1_max', 'x2_max', 'x3_max', 'Gamma', 'Solver_Time', 'Solver_info', 'Method', 'Use_MinLFR' };
    Results = cell2table(cell(0,numel(Cols)), 'VariableNames', Cols);
end

Results = [ Results ; struct2table(s) ];

% save(Results_mat,'Results');
writetable(Results,Results_csv,'Delimiter','|');
clear Results


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
persist.stoplog;

end