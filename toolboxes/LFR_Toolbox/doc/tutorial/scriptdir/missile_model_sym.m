% The continuum of linearized models of a missile is modelled as
% an LFR-object
%
% Symbolic approach.

% Load numerical data
   missiledata

% Define symbolic objects
   syms Al q Ma dp

% Build differential equations
   Cz = z3*Al^3 + z2*Al^2 + z1*( 2 -(1/3)*Ma)*Al + z0*dp;
   Cm = m3*Al^3 + m2*Al^2 + m1*(-7 +(8/3)*Ma)*Al + m0*dp;

   A1 = q+K1*Ma*Cz*(1-Al^2/2);%+Al^4/24);
   A2 = K2*Ma^2*Cm;
   C1 = K3*Ma^2*Cz;

   F = [A1;A2;C1];

% Differentiate for obtaining linearized models
   ABCD = [diff(F,'Al') diff(F,'q') diff(F,'dp')];

% PLug equilibrium surface into ABCD
   dp = -(m3*Al^3 + m2*Al^2 + m1*(-7 +(8/3)*Ma)*Al)/m0;
   ABCD = eval(ABCD);

% Realization
   ABCD = symtreed(ABCD);

% Input/output form
   sys = abcd2lfr(ABCD,2);

% Order reduction after realization
   sys = minlfr(sys,1000*eps);

% Normalisation
   sys = normalizelfr(sys,{'Al','Ma'},[0  2],[0.24  4]);

   sys_sym = sys;
