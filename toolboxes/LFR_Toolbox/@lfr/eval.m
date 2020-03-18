% EVAL         - Evaluates an LFR-object from values in workspace
%-----------------------------------------------------------------
% PURPOSE
% Evaluates an lfr-object from parameter values in the workspace.
%
% SYNOPSIS
% sys2 = eval(sys1)
%
% DESCRIPTION
% - Replaces  the  parameters of the Delta-matrix of an lfr-object
%   by  their  values  as available in the workspace. These values
%   can be numerical values, ss-objects or lfr-objects.
% - Int (i.e., 1/s) can be set to Inf for DC-gain computation.
% - For the LFT: FL([A B;C D],Delta) in which A,B,C or D are given
%   in  the workspace as lfr-objects, the function eval merges the
%   Delta matrices of A, B, C or D with the original one.
% Note that unlike with sym/eval it is not necessary to define all
% parameters in the workspace before using the function eval.
%
% CAUTION:
% In  order  to  avoid conflicts between the internal variables of
% this  function  and  the  parameters  of  the  lfr-object  being
% treated, parameters names shouln't finish with underscore (xx_).
%
% INPUT ARGUMENT
% sys1    lfr-object
%
% OUTPOUT ARGUMENT
% sys2    lfr-object
%
% See also uplft, sym/eval
%-----------------------------------------------------------------
% % EXAMPLE 1: Replacement of a PARAMETER by a function of other
% % parameters and of another one by a linear system.
%
%    lfrs a b c Int
%    sys = [a+b Int/b-c ; 3*Int a^2-b^2];
%
% % Definition of replacements and evaluation
%    a=ss(1,1,1,1); c=b^2;
%    sys = eval(sys);
%
%-----------------------------------------------------------------
% % EXAMPLE 2: Linear system with COEFFICIENTS of A,B,C,D being
% % lfr-objects.
%
%    lfrs a b c Int
%    A = [1+a b*c;2 b]; B = [2*a;5]; C = [1 1]; D = b;
%
% % Prameter dependent linear system
%    blk.names = {'Int'};
%    blk.desc  = [2;2;0;1;1;1;0;0;0;0];
%    sys = lfr(A,B,C,D,blk);
%    size(sys) % 2 states, no uncertainties
%
% % Evaluation of state space matrices
%    sys = eval(sys);
%    size(sys) % 2 states, uncertainties a,b,c
