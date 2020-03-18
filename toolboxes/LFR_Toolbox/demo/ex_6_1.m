% This example illustrates the symbolic approach for modelling
% the continuum of the linearized models of a nonlinear system.
% Here, the considered system is a simplified missile model. First,
% the symbolic expression of the system equations are computed
% (see the manual):

   missiledata
   syms Al q Ma dp

   Cz = z3*Al^3 + z2*Al^2 + z1*( 2 -(1/3)*Ma)*Al + z0*dp;
   Cm = m3*Al^3 + m2*Al^2 + m1*(-7 +(8/3)*Ma)*Al + m0*dp;

   A1 = q+K1*Ma*Cz*(1-Al^2/2);%+Al^4/24);
   A2 = K2*Ma^2*Cm;
   C1 = K3*Ma^2*Cz;

% Now, the state-space matrices of the linearized models are computed
% by deriving the above expressions:

   fg = [A1;A2;C1];
   ABCD = [diff(fg,'Al') diff(fg,'q') diff(fg,'dp')];

% The system matrix  [A B;C D] is realized using the tree
% decomposition (function "symtreed"), then, the input/output
% corresponding form is computed using "abcd2lfr". The result
% is reduced using "minlfr".

   ABCD = symtreed(ABCD);
   sys = abcd2lfr(ABCD,2);
   sys = minlfr(sys,1000*eps);

% The parameter variations must be normalized.

   sys = normalizelfr(sys,{'Al','Ma'},[Al_0-Al_S Ma_0-Ma_S],[Al_0+Al_S Ma_0+Ma_S]);

% Finally, the actuator is added in series at the system input:

   sys = sys*ss(tf([omegaa^2],[1 2*kia*omegaa omegaa^2]));
   sys_6_1 = sys;
