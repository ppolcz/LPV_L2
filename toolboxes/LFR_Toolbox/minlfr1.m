% MINLFR1      - Reduces the order of an LFR
%-----------------------------------------------------------------
% PURPOSE
% Reduces the order of an LFR
%
% SYNOPSIS
% sys2 = minlfr1(sys1[,tol[,iblk[,nocom]]])
%
% DESCRIPTION
% The  technique  behind  is  the  usual  1D - minimal realisation
% method.  This  technique  is applied sequencially to the Dynamic
% block and to uncertainty blocks (given by iblk).
%
% INPUT ARGUMENTS
% sys1   LFR to be simplified.
% tol    Scalar  or  vector of same size as iblk giving tolerances
%        for use of gsminr. Default: tol=1e-7.
% iblk   Vector  numbering the uncertainty blocks to be considered
%        for  order  reduction.  Block number 1 corresponds to the
%        Dynamics.  By  default (iblk = []) all blocks are consid-
%        ered once in natural order.
% nocom  If nocom = 'y' comments are not displayed.
%
% OUTPUT ARGUMENT
% sys2   Reduced LFR.  To check the distance between sys1 and sys2
%        sys2 use distlfr.
%
% See also minlfr, distlfr, lfr/size
%#----------------------------------------------------------------
% % EXAMPLE
%    lfrs a b;
%    sys=1/(a-a+b);
%    size(sys);
%
%    sysr=minlfr1(sys);
%    size(sysr);
%    distlfr(sys,sysr)
