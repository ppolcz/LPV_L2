% CONJ         - Conjugates an LFR-object
%-----------------------------------------------------------------
% PURPOSE
% Computes the conjugate of a lfr-object.
%
% CAUTION
% The function lfr/conj.m  is  not relevant for lfr-objects with
% complex uncertainties or with dynamics.
%
% SYNOPSIS
% sys2 = conj(sys1)
%#----------------------------------------------------------------
% % EXAMPLE
%   sys1 = rlfr(0,2,3,3,5)
%   distlfr(conj(transp(sys1)),sys1')
