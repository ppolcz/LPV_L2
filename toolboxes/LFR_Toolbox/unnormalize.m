% UNNORMALIZE  - Unnormalisation an lfr-object
%-----------------------------------------------------------------
% PURPOSE
% Very  basic provisional unnormalization routine, doesn't work in
% several cases.
%
% SYNOPSIS
% sys_out = normalizelfr(sys);
%
% DESCRIPTION
% sys2 = normalizelfr(sys1) so sys1 = unnormalize(sys2)
% Note  that  if  normalizelfr  has been used with equal lower and
% upper  bounds or radius equal to zero for some parameters, these
% parameters  disappear  from  the block structure so, will not be
% recovered by invoking unnormalize.
%
% INPUT ARGUMENTS
% sys        lfr-object.
%
% OUTPUT ARGUMENT
% sys_out    Returned unnormalized lfr-object.
%
% See also normalizelfr
%#----------------------------------------------------------------
% % EXAMPLE
%    lfrs x y z [1 2 3] [3 4 5] [2.2 3 4.4]
%    lfrs X complex [1+i] [2]
%    M1 = [X+x^2*z (1+y)/(1+z^4) X/(1-z);x^3 1/(2+y)^3 x^3/(2+y)^2];
%    size(M1)
%    M2 = normalizelfr(M1);
%    size(M2)
%    M3 = unnormalize(M2);
%    size(M3)
%    distlfr(M1,M3)
