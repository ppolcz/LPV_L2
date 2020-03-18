% MIN_MAX      - LFR-object: minimum value, maximum value
%-----------------------------------------------------------------
% PURPOSE
% Computes  the  maximum and minimum values of a 1-by-1 lfr-object
% depending on real parameters and without dynamic states.
%
% SYNOPSIS
% [min_val,max_val,min_int,max_int] =
%                 min_max(Mi[,Lambda[,dmin,dmax[,talk[,option]]]])
%
% DESCRIPTION
% Use  first  this  function  with  only  Mi as input argument. If
% there  are  no warnings the result is valid. In case of warnings
% you  must  select  the  second  argument as a gridding of lambda
% (see Method) such that the drawn curve (talk = 2) is less than 1
% on  the  left,  more than 1 in the middle and less than 1 on the
% right.  The  second  argument  can  also  be  used  for a better
% precision, for example try (see Example 2):
% [min_val,max_val] = min_max(Mi);
% [min_val,max_val] = min_max(Mi,[min_val:0.1:max_val]);
%
% METHOD
% This  function  analyses  the well-posedness of inv(Mi - lambda)
% for lambda varying from -inf to +inf.
% - The  first  value  of  lambda  for  which well-posedness fails
%  (i.e.  some mu-measure  becomes larger than one) corresponds to
%  the minimum value of Mi (parameters in the unit ball).
% - The last value of lambda for which well-posedness is satisfied
%  (i.e.  mu  becomes  less  than  one) corresponds to the maximum
%  value of Mi
%
% INPUT ARGUMENTS
% Mi      1-by-1 lfr-object, no dynamics, only real uncertainties.
% Lambda  Vector  of  real values (see Description above), if = []
%         a default value of Lambda is computed.
% dmin    Minimum values of the uncertain parameters ordered as in
%         Mi.blk.names  (excluding the constant block if present).
%         By default, it is the bounds defined in Mi.blk.desc that
%         are considered.
% dmax    Maximum values (as above).
% talk    Integer: 0 no comments, 1: progression of mu computation
%         displayed, 2 (default): the mu-curve is drawn.
% option  String 'default', 'fast' or 'accurate' (very slow). Note
%         that accuracy  is  also  related to the tightness of the
%         gridding defined by the second parameter.
%
% OUTPUT ARGUMENTS
% min_val Minimum  of Mi for parameters variations defined by dmin
%         and dmax.
% max_val Maximum  of Mi for parameters variations defined by dmin
%         and dmax.
% min_int interval for minimal value.
% max_int interval for maximal value.
%
% See also udistlfr, distlfr
%#----------------------------------------------------------------
% % EXAMPLES
%
% % Example 1: a in [-1 1], b in [-1 1]
%    lfrs a b
%    Mi = (10+a*b^2-(1+b))^2/(10+a*b^2);
%    [min_val,max_val] = min_max(Mi,[],[],[],1)
% % Verification of the above bounds by a 3-D plot
%    plotlfr(Mi,{'a',-1,1,10},{'b',-1,1,10})
%
% % Example 2: a in [2,3], b in [3,4], c in [4,5]
% % Bounds computed manually:  19 < Mi <37
%    lfrs a b c [2,3,4] [3,4,5]
%    Mi = 1+a^2+a*b+a*c;
%    [min_val,max_val,min_int,max_int] = min_max(Mi)
% % More precision by zooming around above bounds
%    [min_val,max_val,min_int,max_int] = min_max(Mi,[min_val:0.1:max_val])
