% IMAG         - Imaginary part of LFR-objects
%-----------------------------------------------------------------
% PURPOSE
% Computes the imaginary part of a lfr-object.
%
% CAUTION
% The function lfr/imag.m  is  not relevant for lfr-objects with
% complex uncertainties or with dynamics.
%
% SYNOPSIS
% sys2 = imag(sys1)
%#----------------------------------------------------------------
% % EXAMPLE
%    sysr = rlfr(0,5,7);
%    sysi = rlfr(0,5,7);
%    sys = sysr + j*sysi;
%    distlfr(sysi,imag(sys))
