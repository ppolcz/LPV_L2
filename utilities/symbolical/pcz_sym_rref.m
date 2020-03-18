function [ret] = pcz_sym_rref(A, varargin)
%% pcz_sym_rref
%  
%  File: pcz_sym_rref.m
%  Directory: utilities/symbolical
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2019. June 07. (2019a)
%

args.tolerance = 1e-10;
args = parsepropval(args,varargin{:});

%%

A_rref = rref(A);

% Rounding #1 -- detect zeros
A_rref = A_rref(abs(A_rref) < args.tolerance);

I = find(A_rref);


end