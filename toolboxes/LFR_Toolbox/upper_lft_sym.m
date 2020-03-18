% UPPER_LFT_SYM - Calculates the upper LFT symbolically
%-----------------------------------------------------------------
% PURPOSE
% This  function  can be called to calculate the symbolic transfer
% function  matrix  (upper  LFT).  Only  for systems with repeated
% real or complex blocks. Wrong result or error for full blocks.
%
% SYNOPSIS
% out = upper_lft_sym(sys)
%
% INPUT ARGUMENTS
% sys      dlfr-object
%
% OUTPUT ARGUMENT
% out      symbolic transfer function matrix.
%
%-----------------------------------------------------------------
% %Example
%    dlfrs par1 par2;
%    sys=par1+par2+1/par2;
%    out=upper_lft_sym(sys)
