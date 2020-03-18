% MINLFR       - Reduces the order of an LFR object
%-----------------------------------------------------------------
% PURPOSE
% Calculates   a   minimal,  multi-D  realization  by  eliminating
% unreachable  and  unobservable  subspaces. Method of R. D'Andrea
% and S. Khatri (ACC 1997, pp 3557-3461).
%
% SYNOPSIS
% sysout = minlfr(sys[,tol]);
%
% WARNING
% Choose  by  trials  and  errors the second input argument as its
% default value is often too small.
%
% INPUT ARGUMENTS
% sys    LFR-object.
% tol    Optional tolerance argument, default value see ctrfb.
%
% OUTPUT ARGUMENT
% sysout Returned LFR-object.
%
% See also lfr/size, distlfr, minlfr1
%#----------------------------------------------------------------
% % EXAMPLE
%
%   lfrs a b;
%   sys=[1/(a-a+b) 1/(a-a+b)];
%   size(sys);
%
%   sysr=minlfr(sys);
%   size(sysr);
%   distlfr(sys,sysr)
