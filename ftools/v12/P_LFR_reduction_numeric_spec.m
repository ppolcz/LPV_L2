function [pLFR_reduced] = P_LFR_reduction_numeric_spec(syslfr,varargin)
%% P_LFR_reduction_for_LFR
%  
%  File: P_LFR_reduction_for_LFR.m
%  Directory: 1_PhD_projects/00_my_toolboxes/FinslerTools/v11
%  Author: Peter Polcz (ppolcz@gmail.com) 
% 
%  Created on 2019. March 11.
%  Modified on 2019. March 29. (simplified version)
% 

%% Initialization

args.State = [];
args.Multiplier = []; % TODO: general right LFR multiplier (konnyu kiegeszites)
args.Round = 10;
args = parsepropval(args, varargin{:});

TMP_IbSWJNMuIiKbocfQKqXb = pcz_dispFunctionName;

%% Collect LFR model matrices

pLFR = syslfr;
if isa(pLFR,'lfr')
    pLFR = plfr(syslfr);
end

[A,B,C,D,Delta,bounds] = pLFR.data;

s1 = pLFR.ny;
s2 = pLFR.nu;
m1 = pLFR.m1;
m = s2+m1;

% pcz_display(s1,s2,m1,m);

%% Generate PI and compute its constant left kernel

% [M11,M12] = pcz_split_matrix(eye(m),m,[s2 m1]);
% PI_lfr = plfr(M11,M12,C,D,blk);

PI_lfr = pLFR.generatePI;

if ~isempty(args.State)
    % Nonlinear case: f(x,p) = A(x,p)*x = [F11 F12] * [x;pi_1(x,p)]
    %  where pi(x,p) = [x;pi_1(x,p)] = [I;PI_1(x,p)] * x

    x = args.State;
    one = ones(size(x));
    x_lfr = plfr(one*0,diag(one),one,diag(one*0),diag(x),[-one one]);
    
    PI_lfr = plfr(PI_lfr * x_lfr);
end
    
% LPV case: f(x,p) = A(p)*x
%  where A(p) = [F11 F12] * PI(p) 
%  and PI(p) = [I;PI_1(p)]

rcond_tol = 1e-4;

while true
    %%
    Ker = P_constKer_for_LFR(PI_lfr')';


    Theta = null(Ker);

    % Dimensions
    [~,K] = size(Theta);
    wh_m = rank(Theta);
    wh_m1 = wh_m - s2;

    % pcz_display(s1,s2,wh_m,wh_m1,K);

    %% Generate permutation matrix

    [~,Irows] = rref(Theta');
    Jrows = setdiff(1:m,Irows);
    sigma_b = [Irows Jrows];
    Is_right = pcz_permat(sigma_b)';
    Is1 = Is_right(s2+1:end,s2+1:end);
    Is_left = blkdiag(eye(s1),Is1);

    assert(all(sigma_b(1:s2) == 1:s2), ...
        'The first n elements of sigma_b must be identity permutation (in theory)');

    % pcz_display(Irows,Jrows,sigma_b,Is_right,Is1)

    %% Generate transformation matrix

    V = Theta(Irows,:);
    W = Theta(Jrows,:);

    if rcond(V*V') < rcond_tol
        rcond_tol = rcond_tol / 10;
        pcz_dispFunction('Tolerance value decreased to %g', rcond_tol)
        continue
    end

    Gamma = W*V'/(V*V');

    %{
        args.Round = 10;
    %}

    if args.Round
        Gamma = round(Gamma, args.Round);
    end

    % Model reduction transformation driven by:
    % I_sigma*PI = S_sigma*wh_PI
    % Ss: S_sigma
    Ss = [
        eye(wh_m)
        Gamma  
        ];

    %{    
        pInv_S = [ eye(wh_m) zeros(wh_m,m-wh_m) ]

        PI = sym(PI_lfr)
        wh_PI = pInv_S*Is_right*PI

        Is_right*PI - Ss*wh_PI

        pcz_display(V,W,Gamma,S,pInv_S)
        pcz_display(Is_left,M,Is_right,Ss)

    %}


    %%
    % Model reduction transformationsz

    M = [A B ; C D];

    [A,B,C,D] = pcz_split_matrix(Is_left*M*Is_right'*Ss, [s1 wh_m1],[s2 wh_m1]);
    Delta = pcz_split_matrix(Is1 * Delta * Is1',wh_m1,wh_m1);

    % A,B,C,D,Delta,bounds
    
    pLFR_reduced = plfr(A,B,C,D,Delta,[],pLFR.blk);

    %{
        % TEST:
        simplify(sym(pLFR_reduced) - sym(pLFR))
    %}
    
    %%
    
    break
    
end

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
      
    global SCOPE_DEPTH VERBOSE LATEX_EQNR 
    SCOPE_DEPTH = 0;
    VERBOSE = 1;
    LATEX_EQNR = 0;

    pcz_generateSymStateVector(3,'x')
    A = simplify([ (x1^2 + x2 + 1)/(x3^4 + 3*x1^3 + 4*x3 + 2) + x2 1 1 ]);
    
    warning off
    pLFR = plfr(sym2lfr(A));
    warning on
    
    [pLFR_sym,PI_sym] = sym(pLFR);
    
    pcz_symzero_report(A - pLFR_sym, 'A(x) = pLFR(x)')
    
    %%
    
    pLFR_test = P_LFR_reduction_numeric(pLFR, 'State', x, 'Round', 10);

    [pLFR_test_sym,PI_test_sym] = sym(pLFR_test);

    pcz_symzero_report(A - pLFR_test_sym, 'A(x) = pLFR_test(x)')
    
    %%

    [pLFR_reduced_regi,pLFR_regi] = P_LFR_reduction(A,x);

    [pLFR_regi_sym,PI_regi_sym] = sym(pLFR_regi);

    pcz_symzero_report(A - pLFR_regi_sym, 'A(x) = pLFR_regi(x)')
    
    PI_regi_sym
    PI_test_sym
    
end

function test_reciproc
    %% RECIPROC -- DOES NOT WORK FOR THIS MODEL 
    
    syms a real
    pLFR = plfr(sym2lfr([1;1]/a))
    pLFR_test = P_LFR_reduction_numeric(pLFR, 'Round', 10)

end

function test2
    %% RECIPROC -- DOES NOT WORK FOR THIS MODEL
      
    global SCOPE_DEPTH VERBOSE LATEX_EQNR 
    SCOPE_DEPTH = 0;
    VERBOSE = 1;
    LATEX_EQNR = 0;

    bounds = [
        1 10
        1 10
        ];
    
    [x,x_cell] = pcz_generateLFRStateVector('x',bounds);
    f_fh = @(x1,x2) [
        x1/x2
        x2^2
        ];

    n = numel(x_cell);

    pLFR = plfr(f_fh(x_cell{:}))
    
    [f_sym,Pi_sym] = sym(pLFR)

    pLFR_test = P_LFR_reduction_numeric(pLFR,'Round',10);
   
    [pLFR_test_sym,Pi_test] = sym(pLFR_test)
    
    pcz_symzero_report(f_sym - pLFR_test_sym)
    
end

function test_fermentor
    %% Bioreactor (simple) -- DOES NOT WORK FOR THIS MODEL
        
    global SCOPE_DEPTH VERBOSE LATEX_EQNR 
    SCOPE_DEPTH = 0;
    VERBOSE = 1;
    LATEX_EQNR = 0;
    

    Sf = 10;
    Kr = 1;
    Y = 0.5;
    F = 10;

%     [x,x_cell] = pcz_generateLFRStateVector('x',4);
%     f_fh = @(x1,x2,x3,F) [
%         Kr*x1*x2 - x1/x3*F
%         -Kr/Y*x1*x2 + (Sf-x2)/x3*F
%         F
%         0
%         ];

    [x,x_cell] = pcz_generateLFRStateVector('x',3);
    f_fh = @(x1,x2,x3) [
        Kr*x1*x2 - x1/x3*F
        -Kr/Y*x1*x2 + (Sf-x2)/x3*F
        F
        ];

    n = numel(x_cell);

    pLFR = plfr(f_fh(x_cell{:}))

    [f_sym,Pi_sym] = sym(pLFR)

    f = P_LFR_reduction_numeric_spec(pLFR,'Round',10);

    
end
