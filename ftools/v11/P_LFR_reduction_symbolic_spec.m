function [lfr_reduced] = P_LFR_reduction_symbolic_spec(lfr,x,varargin)
%% P_LFR_reduction
%
%  File: P_LFR_reduction.m
%  Directory: 1_PhD_projects/00_my_toolboxes/FinslerTools/v11
%  Author: Peter Polcz (ppolcz@gmail.com)
%
%  Created on 2018. November 21.
%
%
% Inherited from Script P_LFR_reduction_v7
%
%  File: P_LFR_reduction_v7.m
%  Directory: /1_PhD_projects/00_my_toolboxes/algo_P
%  Author: Peter Polcz (ppolcz@gmail.com)
%
%  Created on 2018. March 21.
%  Modified on 2018. April 06. (Log messages + efficiency)
%
% The second argument should always be the state vector
%
% Example:
%
%   P_LFR_reduction_vN(sym2lfr(A_sym),x)
%
%   P_LFR_reduction_vN(P_lfrdata_vN(sym2lfr(A_sym)),x)
%

%% Initialization

args = struct();
args = parsepropval(args, varargin{:});

TMP_IbSWJNMuIiKbocfQKqXb = pcz_dispFunctionName;

%%

A_sym = [];
if isa(lfr,'sym')

    A_sym = lfr;
    lfr_reduced.LFR_symbolic.f_sym = simplify(A_sym*x);
    lfr_reduced.LFR_symbolic.F_sym = A_sym;
    lfr_reduced.LFR_symbolic.Right_multplier_vector = x;
    try

        warning off
        lfr = sym2lfr(A_sym);
        warning on

        % valamiert a sym2lfr soran elvesziti a valos erteku feltetelezest
        % 2018.09.26. (szeptember 26, szerda), 10:30
        assumeAlso(in(x,'real'));

    catch ex

        pcz_warning(false, 'sym2lfr failed: %s', pcz_getReport(ex));

        vars = symvar(A_sym);

        if isempty(vars)
            lfr_reduced = empty_lfr(lfr_reduced,double(A_sym),x);
            return
        end


        vars_cell = num2cell(vars);

        vars_names = cellfun(@char, vars_cell,'UniformOutput', false);
        names_space = join(vars_names, ' ');

        cmd = sprintf('lfrs %s; vars_lfr = { %s };', names_space{:}, names_space{:});

        pcz_info('Trying to evaluate code: `%s`', cmd);

        vars_lfr = {};
        eval(cmd)

        A_sym_fh = matlabFunction(A_sym, 'vars', vars);
        lfr = A_sym_fh(vars_lfr{:});

        pcz_dispFunction ''
        pcz_info('LFR model build in an other way.')
    end

else
    lfr_reduced = struct;
end

if isfield(lfr,'x') && nargin < 2
    x = lfr.x;
end

if ~isfield(lfr,'A') || ~isfield(lfr,'Delta')
    lfr = P_lfrdata(lfr);
end

tA = lfr.A;
tB = lfr.B;
tC = lfr.C;
tD = lfr.D;
tDelta = lfr.Delta;

% pcz_display(tA,tB,tC,tD,tDelta)

%% Redundant LFR model

tic

[~,n] = size(tA);
[~,p] = size(tB);

% Original model
tF = eye(p) - tDelta*tD;
tG = -tDelta*tC;
tPi = -(tF\tG)*x;
tEq = tA*x + tB*tPi;

if ~isempty(A_sym)
    % pcz_symeq_report(tA - tB/tF*tG,A_sym, 'A + B/(I - Delta*D)*Delta*C == A_sym')
    % pcz_symeq_report(tA*x - tB*(tF\tG)*x,A_sym*x, '(A + B/(I - Delta*D)*Delta*C)x == A_sym*x')
    pcz_symeq_report(tEq,A_sym*x, 'A*x + B*Pi == A_sym*x')
end

pcz_dispFunction
pcz_dispFunction('Redundant LFR model (dim: %d) (symbolic computations). Elapsed time is %.6f', numel(tPi), toc)

%% Canonical decomposition of ${\tt tPi} = \tilde{\pi}_b$

tPib = [ x ; tPi ];
[tTheta,tPi0,q] = P_Pi_canonical_decomp(tPib);

lfr_reduced.LFR_before_reduction.A = tA;
lfr_reduced.LFR_before_reduction.B = tB;
lfr_reduced.LFR_before_reduction.C = tC;
lfr_reduced.LFR_before_reduction.D = tD;
lfr_reduced.LFR_before_reduction.M = [tA tB ; tC tD];
lfr_reduced.LFR_before_reduction.Delta = tDelta;
lfr_reduced.LFR_before_reduction.Pi = tPi;

lfr_reduced.LFR_before_reduction.Pi_candecomp.Theta = tTheta;
lfr_reduced.LFR_before_reduction.Pi_candecomp.Pi0 = tPi0;
lfr_reduced.LFR_before_reduction.Pi_candecomp.q = q;

%%
% Dimensions
[m,K] = size(tTheta);
nk = rank(tTheta);
k = nk - n;
% pcz_display(n,p,k,K);

%% Permutation of variables of LFR
% Permutation matrix ${\tt Ir} = I_\varrho$.
[~,Icols] = rref(tTheta);
Jcols = setdiff(1:K,Icols);
rho = [Icols Jcols];
Ir = pcz_permat(rho)';

%%
% Permutation matrix ${\tt Isb} = I_{\sigma_b}$ and ${\tt Is} = I_\sigma$.
[~,Irows] = rref(tTheta');
Jrows = setdiff(1:m,Irows);
sigma_b = [Irows Jrows];
Isb = pcz_permat(sigma_b)';
Is = Isb(n+1:end,n+1:end);

% pcz_display(Irows,Jrows,sigma_b,Isb,Is)

%%
% Permuted coefficient matrix
pTheta = tTheta(sigma_b,rho);

% pcz_display(pTheta)

assert(all(sigma_b(1:n) == 1:n), ...
    'The first n elements of sigma_b must be identity permutation (in theory)');
assert(rank(pTheta(1:nk,1:nk)) == nk, ...
    'The upper left (n+k)x(n+k) submatrix of Theta should be full-rank');

%%
% Auxiliary matrices
V = pTheta(1:nk,:);
W = pTheta(nk+1:end,:);
Gamma = round(W*V'/(V*V'),10);

Gamma1 = Gamma(:,1:n);
Gamma2 = Gamma(:,n+1:n+k);

S = Isb' * [
    eye(n+k)
    Gamma
    ];

% pcz_display(V,W,Gamma)

%%
% Model reduction transformationsz

% Isb
% [tA tB ; tC tD]
% Isb'
% S
% tDelta

warning 'FIGYELEM itt meg nincs megoldva a permutacios matrix'

M = pcz_split_matrix(Isb*[tA tB ; tC tD]*Isb'*S,n+k,n+k);
[A,B,C,D] = pcz_split_matrix(M, [n k],[n k]);

Delta = pcz_split_matrix(Is * tDelta * Is',k,k);

Pi0 = Ir * tPi0;

%%

I = eye(k);

F = I - Delta*D;
G = -Delta*C;
Pi = pcz_split_matrix(Is*tPi,k,1);
Pib = [ x ; Pi ];
Eq = A*x + B*Pi;

% pcz_2basews

pcz_symeq(S*Pib,tPib, sprintf('Pi [dim: %d] = S * (reduced Pi [dim: %d]).', numel(tPib), numel(Pib)));

pcz_symeq_report(tEq,Eq,1e-5,'Equation of the reduced model should be the same')

LFR_reduced = struct;
lfr_reduced.LFR_reduced = pcz_struct_append(LFR_reduced,...
    A,B,C,D,M,Delta,Pi,S,Gamma);

%%

pcz_dispFunction
pcz_dispFunctionEnd(TMP_IbSWJNMuIiKbocfQKqXb);

end


function lfr = empty_lfr(lfr,A,x)
    lfr.A = A;
    lfr.B = zeros(size(A,1),0);
    lfr.C = zeros(0,size(A,2));
    lfr.D = zeros(0,0);
    lfr.Delta = zeros(0,0);
    lfr.Pi = zeros(0,1);
    lfr.Pib = x;
end



function test1
%%

syms x real

A = (x^2 + x + 1)/(x^4 + 3*x^3 + 4*x + 2) + x;

plfr_test = P_LFR_reduction_symbolic(A,x)
pLFR = P_LFR_reduction(A,x)

simplify(pLFR.LFR_reduced.Pi - plfr_test.LFR_reduced.Pi)

pLFR.LFR_reduced.M
plfr_test.LFR_reduced.M

pLFR.LFR_reduced.D
plfr_test.LFR_reduced.D

end
