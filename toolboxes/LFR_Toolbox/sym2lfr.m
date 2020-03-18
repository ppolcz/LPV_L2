% SYM2LFR      - From symbolic expression to LFR-object
%-----------------------------------------------------------------
% PURPOSE
% Systematic  object-oriented  realization from a symbolic expres-
% sion.  Preserves  the  parenthesis  structure  of  the  symbolic
% expression  to  be  treated,  normalizes all coefficients of the
% realization (if tol>=0) and remove negligible monomials (tol>0).
%
% SYNOPSIS
% sys = sym2lfr(symbex[,tol])
%
% INPUT ARGUMENTS
% symbex     Symbolic  expression.
% tol        tolerance  argument (default tol = 0).  After scaling
%            coefficients  in  each  parenthesis,  nomomials  with
%            coefficients  of  magnitude  less  than  tol  can  be
%            removed.  Note  that tol > 0 (enabling approximation)
%            is relevant ONLY if all symbolic variables apprearing
%            in symbex belong to the same interval (e.g. [-1 +1]).
%            If tol = -1, there is no factorization/normalization/
%            approximation (equivalent to sys = eval(symbex) after
%            parameters are defined in the worksapce using lfrs).
%
% OUTPUT ARGUMENT
% sys        lfr-object correponding to symbex
%
% See also str, symtreed
%#----------------------------------------------------------------
% % EXAMPLE
%   syms a b c
%   M = [1e-40*a*b*(2*1e+10+3*1e+10*a+4*1e+10*c*b)^2*(1+6*1e+10*c)^2 a^3*b;1 1/(b+c)];
%
% % Transformation to lfr-object with coefficient normalization
%   M1 = sym2lfr(M);
%
% % Transformation to lfr-object without normalization
%   lfrs a b c
%   M0 = eval(M); % <=> M0 = sym2lfr(M,-1);
%
% % Check internal values
%   [norm(M1.a) norm(M0.a);norm(M1.b) norm(M0.b);norm(M1.c) norm(M0.c)]
%   distlfr(M1,M0)
