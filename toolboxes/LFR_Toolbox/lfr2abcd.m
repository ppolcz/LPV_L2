% LFR2ABCD     - From input-output to state-space model
%-----------------------------------------------------------------
% PURPOSE
% Converse of abcd2lfr: From input-output to state-space model.
%
% SYNOPSIS
% [abcd,netats] = lfr2abcd(sys);
%
% INPUT ARGUMENTS
% sys     Corresponding input/output lfr-object.
%
% OUTPUT ARGUMENT
% abcd     lfr-object representing the matrix [a b;c d] of a LTI
%          system (dx/dt = ax + bu, y = cx + du).
% nstates  Number of states
%
% See also  abcd2lfr
%#----------------------------------------------------------------
% % EXAMPLE
% % From [A B;C D] to I/O lfr-object
%   lfrs x4
%   abcd0 = (1/x4)*rlfr(0,8,7,2,4,1);
%
% % From input/output to state-space
%   sys = abcd2lfr(abcd0,5);
%
% % Converse
%   [abcd,netats] = lfr2abcd(sys);
%
% % Comparison
%   distlfr(abcd0,abcd)
