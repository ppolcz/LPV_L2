% This example illustrates the realization of LFR-objects subject
% to complex scalar repeated uncertainties. For example:
%
% ((a + 1/s)/(1 + b/s))*((1+2*(d1+d2)+d1*d2)/(1-2*(d1+d2)+d1*d2))
%
% in which a and b are real, d1 and d2 are complex scalars.

   lfrs  Int a b
   lfrs  d1 d2 'complex'

   sys = (1+2*(d1+d2)+d1*d2)/(1-2*(d1+d2)+d1*d2);
   sys = (a+Int)/(1+b*Int)*sys;

   size(sys)

