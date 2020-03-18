% LF2LFR       - LFR of left factorized representations
%-----------------------------------------------------------------
% PURPOSE
% This function computes the lfr-object SYS with transfer matrix G
% from  the factors N and D of a left  factorization G = D^{-1}*N.
% The factors are specified through the compound system ND = [N D].
%
% SYNOPSIS
% sys=lf2lfr(ndh);
%
% INPUT ARGUMENTS
% ND     lfr-object  containing  the  factors  N  and  D of a left
%        factorization  G = D^{-1}*N; the factors are specified by
%        the column concatenated representation ND = [N D].
%
% OUTPUT ARGUMENT
% sys    The resulting lfr-object with transfer matrix G.
%
% See also rf2lfr
%#----------------------------------------------------------------
% % EXAMPLE
%    lfrs a b c
%    ND = [a^2*[1 1 b];[c 1 1]/b^2];
%    sys1 = lf2lfr(ND);                 % Preserves the size of ND
%    sys2 = inv(ND(:,2:3)) * ND(:,1);   % Size multiplied by 2
