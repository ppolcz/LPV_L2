% LFRTOL       - For systematic order reduction
%-----------------------------------------------------------------
% PURPOSE
% If  the  input  argument of this function is non-zero, the order
% reduction  function  minlfr is invoked systematically each times
% an  addition, multiplication or concatenation is applied to lfr-
% objects.
%
% SYNOPSIS
% lfrtol            <- initializes a global tolerance if necessary
% lfrtol(tol)       <- defines tolerance for nD-reduction (minlfr)
% lfrtol('info')    <- for information
% lfrtol('default') <- default: no systematic order reduction
%
% INPUT ARGUMENTS
% tol     Positive value used as tolerance argument of minlfr.
%
% See also minlfr, @lfr
%#----------------------------------------------------------------
% % EXAMPLES
%    clear global
%    lfrtol('i')
%    lfrtol
%    lfrtol('i')
%    lfrtol(0.0001)
%    lfrtol('i')
