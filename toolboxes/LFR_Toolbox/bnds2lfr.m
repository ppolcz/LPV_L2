% BNDS2LFR     - From matrix bounds to LFR-object
%-----------------------------------------------------------------
% PURPOSE
% Transforms  a  matrix  with  bounds relative to its entries to a
% lfr-object.
%
% SYNOPSIS
% sys = bnds2lfr(names,matmax)
% sys = bnds2lfr(names,matmin,matmax);
% sys = bnds2lfr(names,matnom,percent,'p');
%
% DESCRIPTION
% Let us consider an example. The matrix [a11 a12;a21 a22] is such
% that   a11m < a11 < a11M,  a21m < a21 < a21M,  a21m < a21 < a21M
% and a22 fixed.  For  transformation  to a lfr-object, the input
% arguments must be:
%    matmin = [a11m a12m;a21m a22];
%    matmax = [a11M a12M;a21M a22];
%
% INPUT ARGUMENTS
% name      String.  If  for  example  name = 'Y_', the  uncertain
%           parameters will be named  'Y_11','Y_12' and so on.
% matmin    Real matrix giving elementwise the lower bounds.
% matmax    Similar for upper bounds. If matmin is omitted, matmin
%           = - matmax (so, all entries of matmax must be > 0).
% matnom    Nominal real values (used with option 'percent').
% percent   Same size as the considered uncertain matrix. It gives
%           elementwise the uncertainty percentages. This argument
%           can also be 1_by_1.
% option    Should be 'p' if percentages are considered.
%
% OUTPUT ARGUMENT
% sys       Resulting lfr-object.
%
% See also abcd2lfr, data2lfr
%#----------------------------------------------------------------
% % EXAMPLE 1 (interval of variations)
%    matmin = [0.5 2.0 2.9; 3.9 4.7 5.5];
%    matmax = [1.5 2.0 3.1; 4.1 5.3 5.5];
%    A = bnds2lfr('A_',matmin,matmax); size(A)
%    B = bnds2lfr('B_',matmin); size(B)
%
% % EXAMPLE 2 (pencent)
%    matnom = [0.5 2.0 2.9; 3.9 4.7 5.5];
%    percent = 10;
%    A = bnds2lfr('A',matnom,percent,'p'); size(A)
%
% % EXAMPLE 3 (pencent / matrix)
%    matnom = [0.5 2.0 2.9; 3.9 4.7 5.5];
%    percent = [0 0 10;2 2 0];
%    A = bnds2lfr('A_',matnom,percent,'p'); size(A)
