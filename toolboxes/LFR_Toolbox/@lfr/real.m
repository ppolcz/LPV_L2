% REAL         - Real part of LFR-objects
%-----------------------------------------------------------------
% PURPOSE
% Computes the real part of a lfr-object.
%
% CAUTION
% The function lfr/real.m  is  not relevant for lfr-objects with
% complex uncertainties or with dynamics.
%
% SYNOPSIS
% sys2 = real(sys1)
%#----------------------------------------------------------------
% % EXAMPLE
%    sysr = rlfr(0,5,7);
%    sysi = rlfr(0,5,7);
%    sys = sysr + j*sysi;
%    distlfr(sysr,real(sys))
