% DISTLFR      - Distance between two LFR-objects
%-----------------------------------------------------------------
% PURPOSE
% Computes a "distance" between two lfr-objects
%
% SYNOPSIS
% [dist1,dist2,dist3,dist4] = distlfr(lfr1,lfr2[,nrand[,frequ]])
%
% DESCRIPTION
% This  function computes a distance measure by evaluating (lfr1 -
% lfr2) at random values.
% The  uncertain  parameters are chosen at random on the unit ball
% and  the  frequency w (or z for discrete-time systems) is chosen
% in [0 1] unless the fourth input argument  specifies alternative
% bounds. A  specific treatment is considered  for  real  repeated
% uncertainties: as a priority, depending on the  chosen number of
% random  trials (nrand), the  difference  between  both  compared
% objects is considered at all or part of the vertices of the unit
% hypercube. For  real repeated uncertainties with min/max-bounds,
% the random-values are chosen within the interval [min,max].
% The  returned  ``distances''  are  relative  or  absolute. These
% distances  are  also  available entry per entry. Relative values
% should be understood as follows:  If the  considered lfr-objects
% are normalized so that the uncertain parameters vary in [-1,+1],
% a returned value of  dist1 = 0.05  means that the worst transfer
% difference is about 5 percent.
%
% CAUTION
% If one entry of lfr1 or lfr2 is indentically equal to zero, the
% relative  distance (= dist1) is equal to 1. In this case, it is
% more relevant to consider the absolute distance (= dist2). More
% generally,  consider  the  absolute  distance  any  times  some
% returned relative value is exactly equal to 1.
%
% INPUT ARGUMENTS
% lfr1,lfr2 Objects to be compared.
% nrand     Number of random trials (default: 200).
% frequ     1-by-2  matrix  giving frequ. bounds (default: [0 1]).
%
% OUTPUT ARGUMENTS
% dist1     Worst  relative "distance" (maximum value in dist3).
% dist2     Worst  absolute "distance" (maximum value in dist4).
% dist3     Worst  relative "distance" entry per entry.
% dist4     Worst  absolute "distance" entry per entry.
