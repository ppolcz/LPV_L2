% REMALGEQ     - removes redunant algebraic equations
% ----------------------------------------------------------------
% PURPOSE
% Constant blocks are introduced in order to make feasible
% inversion of lfr-objects with non-invertible D-matrix. When lfr-
% objects  are  evaluated  these blocks might become useless. This
% function reduces the constant block to its minimum size.
%
% SYNOPSIS
% sys2 = remalgequ(sys1)
%
% REMARK
% This  function  can be ignored, it is invoked by other functions
% of the toolbox (minlfr, uplft, eval)
% ----------------------------------------------------------------
