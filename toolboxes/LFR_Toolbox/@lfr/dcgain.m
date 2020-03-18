% DCGAIN       - Steady-state gain of an LFR-object
%-----------------------------------------------------------------
% PURPOSE
% Computes the steady-state gain of a dynamic lfr-object
%
% SYNOPSIS
% g = dcgain(sys)
%#----------------------------------------------------------------
% % EXAMPLES
%   lfrs Delay Int a b c
%
%   sys1 = (1+a+b*Int)/(1+c+Int)
%   g1 = dcgain(sys1)
%
%   sys2 = (1/a+b*Delay+c*Delay^3)/(1+Delay+a*b*Delay^2);
%   g2 = dcgain(sys2)
