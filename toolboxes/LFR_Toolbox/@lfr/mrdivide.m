% MRDIVIDE     - LFR-object division
%--------------------------------------------------------------------
% PURPOSE
% Computes  the  inverse  of an lfr if invoked by 1/sys, or computes
% a division if invoked by sys1/sys2 ( = sys1*(1/sys2) ).
%
% SYNOPSIS
% sys = 1/sys1;
% sys = sys1/sys2;
%
% DESCRIPTION
% The  inverted lfr must have the same number of inputs and outputs.
% If sys2.d is not invertible, a generalized inverse (descriptor-LFT)
% is calculated. In the division case (sys=sys1/sys2),  the number of
% inputs of sys1  must  be equal to the number of inputs of sys2.
%
% See also lfr/minus, lfr/uminus, lfr/mtimes, lfr/plus
%--------------------------------------------------------------------
