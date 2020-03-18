% UDISTLFR     - Upper-bound of LFR-objects distance
%-----------------------------------------------------------------
% PURPOSE
% Computes an upper-bound of the distance between two lfr-objects.
%
% SYNOPSIS
% [dist1,dist2,mindiff,maxdiff]=udistlfr(lfr1,lfr2[,talk]]);
%
% DESCRIPTION
% This function computes the distance between two real lfr-objects.
% It treats separately each entry of the considered lfr-object. It
% is in fact an  UPPER bound of this distance that is computed. In
% order  to  have a LOWER bound, use the function distlfr (usually
% much  faster for lfr-objects having small number of real uncert-
% ain arameters).
% This  upper bound computation problem is viewed as an artificial
% mu-analysis problem.
% This  function  is useful for evaluating precisely approximation
% errors  (e.g.,  after use of reduclfr). Approximation errors can
% be transformed to lfrs using bnds2lfr (see the example below).
%
% INPUT ARGUMENTS
% lfr1,lfr2  lfr-objects  with  only  real  uncertainties and real
%            ceofficients.
% talk       Integer = 0,1 or 2. No comments: 0, default: 1.
%
% OUTPUT ARGUMENTS
% dist1      Absolute ``distance'' (maximum value in dist2).
% dist2      Absolute ``distance'' entry per entry.
% mindiff    lfr1 = lfr2 + error s.t.:
% maxdiff    mindiff(i,j) < error(i,j) < maxdiff(i,j)
%
% See also min_max, distlfr, bnds2lfr
%#----------------------------------------------------------------
% % EXAMPLE
%     lfrs a b c
%     sys1 = [3+0.001*a^5-b*c a^4*c;a*b*c^3+2 a^2-0.001*b*c^3];
%
% % Example 1: Distance of sys1 from zero
%     sys2 = zeros(2,2);
%     % Lower bound
%     [x,dist_low] = distlfr(sys1,sys2);
%     % Upper bound
%     [dist_up,dist2,mindiff,maxdiff] = udistlfr(sys1,sys2);
%     dist_low
%     dist_up
%
% % Example 2: use of udistlfr and bnds2lfr for approximation error
% % modelling. First, approximation of sys1:
%     sys2 = reduclfr(sys1,0.01,'a');
%
% % Bounds of approximation error
%     [distu,dist2,mindiff,maxdiff] = udistlfr(sys1,sys2);
%
% % Replacement of approximation error by additional real parameters
%     syserr = bnds2lfr('approxerr_',mindiff,maxdiff);
%     sys3 = sys2 + syserr;
%     size(sys3)
