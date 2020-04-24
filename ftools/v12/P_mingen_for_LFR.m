function [S,syslfr_min,iS,Ker] = P_mingen_for_LFR(syslfr,varargin)
%% P_mingen_for_LFR
%  
%  File: P_mingen_for_LFR.m
%  Directory: 1_PhD_projects/00_my_toolboxes/FinslerTools/v11
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2019. March 11.
%

%% Initialization

args.State = [];
args.Round = 10;
args.proj = [];
args.lims = [];
args = parsepropval(args, varargin{:});

TMP_EOLHIWBuIwgwCuuNBzen = pcz_dispFunctionName;

%%

PI = syslfr;
if isa(PI,'lfr')
    PI = plfr(syslfr);
end
Pi = PI;

nu = PI.nu;
m = PI.ny;

if ~isempty(args.State) && isa(args.State,'sym')
    error 'TODO: EZ MEG NEM JO'
    % Nonlinear case: f(x,p) = A(x,p)*x = [F11 F12] * [x;pi_1(x,p)]
    %  where pi(x,p) = [x;pi_1(x,p)] = [I;PI_1(x,p)] * x

    x = args.State; %#ok<UNRCH>
    one = ones(size(x));
    x_lfr = plfr(one*0,diag(one),one,diag(one*0),diag(x),[-one one]);
  
    PI = Pi;
    Pi = plfr(Pi * x_lfr);
end

if ~isempty(args.State) && isa(args.State,'lfr')
    % Nonlinear case: f(x,p) = A(x,p)*x = [F11 F12] * [x;pi_1(x,p)]
    %  where pi(x,p) = [x;pi_1(x,p)] = [I;PI_1(x,p)] * x

    x = args.State;
  
    PI = Pi;
    Pi = plfr(Pi * x);
end

%%

rcond_tol = 1e-7;

%%

Ker = P_constKer_for_GSS(Pi','proj',args.proj,'lims',args.lims)';

Theta = null(Ker);

%% Find V and W, such that V is invertible

% Tolerance to rref (for rank decision).
rank_tol = max(size(Theta)) * eps(1); % norm(Theta) = 1
while true

    [~,Irows] = rref(Theta', rank_tol);
    V = Theta(Irows,:);

    if rcond(V*V') < rcond_tol
        % sigma = svd(Theta)';
        % pcz_dispFunction_num2str(sigma)
        
        rank_tol = rank_tol * 10;
        pcz_dispFunction('Rank decision tolerance for `rref'' increased to %g', rank_tol)
        continue
    end
    
    break
end


Jrows = setdiff(1:m,Irows);
W = Theta(Jrows,:);

sigma = [Irows Jrows];
Is = pcz_permat(sigma)';

% pcz_display(Theta,V,W,Irows,Jrows,sigma)

%% Generate transformation matrix

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
    eye(size(Gamma,2))
    Gamma  
    ];

S = Is'*Ss;

iS = [ eye(size(Gamma,2)) 0*Gamma' ] * Is;

syslfr_min = plfr(iS * PI);
% syslfr_min = syslfr_min.set_vars(PI.subsvars)

%{    

    sym(Pi) - S*sym(syslfr_min)

    pcz_display(V,W,Gamma,S,pInv_S)
    pcz_display(Is_left,M,Is_right,Ss)

%}


pcz_dispFunctionEnd(TMP_EOLHIWBuIwgwCuuNBzen);

%% Igy lehetne sokkal rovidebben
% 
% Ker = P_constKer_for_LFR(syslfr')';
% 
% S = null(Ker);
% 
% iS = pinv(S);
% 
% syslfr_min = minlfr(iS*syslfr);

end