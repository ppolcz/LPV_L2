% DIFF         - Differentiates an LFR-object
%-----------------------------------------------------------------
% PURPOSE
% Differentiate  n  times  an  LFR-object with respect to a single
% parameter
%
% SYNOPSIS
% df = diff(f[,name[,n]])
%
% INPUT ARGUMNENTS
% f           Lfr-object to be differentiated
% name        String :  name  of the parameter for differentiation
%             (cannot  be  Int,  Delay or ConstBlock). If there is
%             only one parameter, this argument can be omitted.
% n           Positive integer : differentiation order. By default
%             n = 1.
%#----------------------------------------------------------------
% % EXAMPLES
% % Note that bounds are defined for well-posedness in distlfr
%    lfrs a b c [1 1 1] [2 2 2]
%    F1 = (1+a+b*c)/(a^2+c*b);
%    dF1 = diff(F1,'a',3);
%
% % Symbolic approach
%    syms a b c
%    F2 = (1+a+b*c)/(a^2+c*b);
%    dF2 = diff(F2,'a',3);
% % From symbolic to LFR plus variation bounds
%    dF2 = sym2lfr(dF2);
%    set(dF2,{'a','b','c'},{[1 2],[1 2],[1 2]},{'minmax','minmax','minmax'});
%
% % Comparison
%    distlfr(dF1,dF2)
