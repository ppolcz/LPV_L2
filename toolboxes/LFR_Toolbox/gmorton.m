% GMORTON      - LFR-object realization
%-----------------------------------------------------------------
% PURPOSE
% LFR-object  realization  by  a generalized Norton approach. This
% function  computes  an  lfr-object  (from  inputs to outputs) in
% which the matrix [A B;C D]  is  given as an expansion.
%
% SYNOPSIS
% [sys,abcd] = gmorton(coeff,expr[,tol])
%
% DESCRIPTION
% Considered dependency:
% - dynamic system
% | A  B |   | A1  B1 |             | A2  B2 |
% | C  D | = | C1  D1 | * expr(1) + | C2  D2 | * expr(2) +...
%
% - non-dynamic matrix
% M = delta1 M1 * expr(1) + M2 * expr(2) + ...
%
% Example:  M = M1 + a*M2 + a*b*M3
%
% => syms a b (or lfrs a b)
%    coeff  = {M1 M2 M3};   (coefficients of the expansion)
%    expr   = [1  a  a*b]   (expression of the expansion)
%
% INPUT ARGUMENTS
% coeff    Cell containing the coefficients (matrices or ss-
%          objects) of the above expansion
% expr     Symbolic or LFR vector.
% tol      Tolerance for SVD decomposition.  Choose this parameter
%          by  trials  and  errors  for  LFT approximation. Defaut
%          value = eps*1e+5.
%
% OUTPUT ARGUMENT
% sys     Resulting LFR-object.
%
% See also abcd2lfr, bnds2lfr, syl2lfr, str
%#----------------------------------------------------------------
% % EXAMPLE
%    sys1 = rss(3,2,3); sys2 = rss(3,2,3); sys3 = rss(3,2,3);
%    M1 = [sys1.a sys1.b;sys1.c sys1.d];
%    M2 = [sys2.a sys2.b;sys2.c sys2.d];
%    M3 = [sys3.a sys3.b;sys3.c sys3.d];
%
% % Use of gmorton with ss-objects in 'coeff'
%    syms a b
%    expr = [1 a b];
%    coeff{1} = sys1; coeff{2} = sys2; coeff{3} = sys3;
%    lfr1 = gmorton(coeff,expr);
%
% % Alternative use of gmorton with matrices in 'coeff'
%    coeff{1} = M1; coeff{2} = M2; coeff{3} = M3;
%    lfri = gmorton(coeff,expr);
%    lfr2 = abcd2lfr(lfri,3);
%
% % Comparison of results
%    distlfr(lfr1,lfr2)
