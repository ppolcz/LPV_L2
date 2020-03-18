% ACTUALVAL    - Actual parameter value after normalization
%-----------------------------------------------------------------
% PURPOSE
% Normalization  using  the  function  normalizelfr applies a non-
% linear  transformation if the nominal value is not in the middle
% of  the range of variation. This function inverts this nonlinear
% transformation for recovering actual parameter values from their
% nomalized values.
%
% SYNOPSIS
% [val] = actualval(sys1,par_name,par_val[,flag])
%
% INPUT ARGUMENTS
% sys        lfr-object  which  has  already been normalized using
%            the function normalizelfr.
% par_name   Cell  of  strings:  names of the parameters for which
%            which we look for actual value.
% par_val    Vector of the normalized numerical values.
% flag       If flag == 'r', actual values are normalized.
%
% OUTPUT ARGUMENT
% act_val    Actual values corresponding to par_val (or normalized
%            values if flag = 'r').
%
% See also normalizelfr
%#----------------------------------------------------------------
% % EXAMPLE
%    lfrs x y z [3 4.5 1] [6 6 5] [4 5 2]
%    lfrs X complex [1+i] [3]
%    M = [x*y X*z^2];
%    size(M)
%    Mn = normalizelfr(M);
%    size(Mn)
%
% % From normalized to actual
%    % maximum values -> [6 6 5]
%    actualval(Mn,{'x','y','z'},[1 1 1])
%    % minimum values -> [3 4.5 1]
%    actualval(Mn,{'x','y','z'},[-1 -1 -1])
%    % nominal values -> [4 5 2]
%    actualval(Mn,{'x','y','z'},[0 0 0])
%    % 1+0i -> 4+i
%    actualval(Mn,{'X'},1)
%
% % normalized to actual to normalized
%    xyz_norm = [1.1 2.2 3.3 4.4i];
%    xyz_act  = actualval(Mn,{'x','y','z','X'},xyz_norm)
%    actualval(Mn,{'x','y','z','X'},xyz_act,'r')     % = xyz_norm
