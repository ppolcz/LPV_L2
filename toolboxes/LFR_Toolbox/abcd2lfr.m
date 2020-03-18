% ABCD2LFR     - LFR of system matrix -> input/output lfr-object
%-----------------------------------------------------------------
% PURPOSE
% Let  abcd  = [a b;c d] (lfr-object) with outputs [d/dt x; y] and
% inputs  [x;u].  The function abcd2lfr transforms abcd to an lfr-
% object  with output y and input u, i.e. the integrator block I/s
% or the delay block I/z are included in the matrix Delta.
%
% SYNOPSIS
% sys = abcd2lfr(abcd,nstates,cont);
%
% INPUT ARGUMENTS
% abcd     LFR-object  representing  the matrix [a b;c d] of a LTI
%          system
% nstates  Number of states
% cont     Continuous-time system(cont=1)(default) or
%          Discrete-time system(cont=0)
%
% OUTPUT ARGUMENT
% sys     Corresponding input/output lfr-object.
%
% See also  lf2lfr, rf2lfr, lfr2abcd
%#----------------------------------------------------------------
