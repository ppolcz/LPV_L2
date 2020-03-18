% RLFR         - Stable random LFR-object with given size
%-----------------------------------------------------------------
% PURPOSE
% Generation  of  a random lfr-object with given number of dynamic
% states (n),  outputs (p), inputs (m),  repetition  of  uncertain
% parameters (n1,n2,...).
%
% SYNOPSIS
% sys = rlfr(n,p,m[,n1[,n2,...]][,prefix])
% sys = rlfr(p,m,blk)
% sys = rlfr(p,m)
%
% INPUT ARGUMENTS
% p       Number of outputs
% m       Number of inputs
% n       Number of states (default: n = 5)
% ni      Size of ith real/scalar block (default: ni = 2)
% prefix  For names of uncertainties (default: prefix = 'x')
% blk     Structured array (see help lfr)
%
% OUTPUT ARGUMENT
% sys     Random lfr-object
%
% See also lfr, lfrs
%#----------------------------------------------------------------
% % EXAMPLE
%    sys1= rlfr(5,2,3,1,0,2,0);
%    size(sys1);
%
%    sys2= rlfr(0,2,3,1,2,'Y');
%    size(sys2);
%
%    blk.names = {'x','Y'};
%    blk.desc  = [[3;3;1;1;1;1;1;2;-1;1;0],[3;5;0;0;1;1;1;2;-2;2;0]];
%    sys3 = rlfr(2,3,blk);
%    size(sys3);
