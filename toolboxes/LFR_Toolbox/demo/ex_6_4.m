% Example ex_6_4 illustrates how to take into account the
% dependencies on the equilibrium surface.

% First, run ex_6_2.m

   ex_6_2

% The first step, in order to take equilibrium dependency, consists
% of writing dp as an LFR-object (see the manual):

   lfrs Al Ma [0.0 2.0] [0.349 4.0]
   lfrs dp
   dp = -(m3*Al^3 + m2*Al^2 + m1*(-7 +(8/3)*Ma)*Al)/m0;

% It remains to replace  dp in the Delta-block by invoking
% "lfr/eval" )

   sys = eval(sys_nn);
   sys = minlfr(sys,1000*eps);

% Normalization (bounds defined when lfrs was invoked)

   sys = normalizelfr(sys);
   sys_6_4 = sys;
   size(sys_6_4)

% This result can be compared to the one of ex_6_3.m
% distlfr(sys_6_3,sys_6_4)
