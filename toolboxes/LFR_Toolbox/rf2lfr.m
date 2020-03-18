% RF2LFR       - LFR of right factorized representations
%-----------------------------------------------------------------
% PURPOSE
% This function computes the lfr-object SYS with transfer matrix G
% from  the factors N and D of a right factorization G = N*D^{-1}.
% The factors are specified through the compound system ND = [N;D].
%
% SYNOPSIS
% sys=rf2lfr(nd);
%
% INPUT ARGUMENTS
% ND     lfr-object  containing  the  factors  N  and D of a right
%        factorization  G = N*D^{-1}; the factors are specified by
%        the column concatenated representation ND = [N ; D].
%
% OUTPUT ARGUMENT
% sys    The resulting lfr-object with transfer matrix G.
%
% See also lf2lfr, abcd2lfr
%#----------------------------------------------------------------
% % EXAMPLE
%    lfrs a b c
%    ND = [a^2*[1;1;b] [c;1;1]/b^2];
%    sys1 = rf2lfr(ND);                 % Preserves the size of ND
%    sys2 = ND(1,:) * inv(ND(2:3,:));   % Size multiplied by 2
