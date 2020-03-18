% REDUCLFR     - Approximation of LFR-objects
%-----------------------------------------------------------------
% PURPOSE
% Handle by a dichotomy search, the size of the tolerance argument
% of  minlfr,  minlfr1...  for a given maximum approximation error
% bound.
%
% SYNOPSIS
% sysout = reduclfr(sysin,errbnd,errtyp[,method[,frequ]]);
%
% INPUT ARGUMENTS
% sysin     Objects to be approximated.
% errbnd    Scalar = approximation error bound.
% errtyp    String 'a' for absolute error, 'r' for relative error.
% method    String = 'minlfr' (default) or 'minlfr1'
% frequ     1-by-2  vector  for  the  frequency range in which the
%           distance is computed ([0 1] by default)
%
% OUTPUT ARGUMENTS
% sysout    Approximation of sysin.
%
% See also minlfr, minlfr1, redlfr1, min_max, distlfr
%#----------------------------------------------------------------
% % EXAMPLE
% % System having a small part
%     lfrs a b c d
%     small = 0.001;
%     sysin = [a 2;c d]*b + small*a^5*[b^3 c;d 1];
%
% % Order reduction by approximation
%     sysout = reduclfr(sysin,0.01,'a');
%     size(sysin);
%     size(sysout);
%     distlfr(sysin,sysout)
