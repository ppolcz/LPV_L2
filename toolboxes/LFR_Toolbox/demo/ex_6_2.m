% Example ex_6_1.m illustrates modelling of a nonlinear system using
% the symbolic approach. Here is illustrated the object-oriented
% counter part. The main function is "difflfr". See the manual
% for explanations.

   missiledata
   lfrs Al Ma [0.0 2.0] [0.349 4.0]  % Bounds specified
   lfrs q dp

   Cz = z3*Al^3 + z2*Al^2 + z1*( 2 -(1/3)*Ma)*Al + z0*dp;
   Cm = m3*Al^3 + m2*Al^2 + m1*(-7 +(8/3)*Ma)*Al + m0*dp;

   A1 = q+K1*Ma*Cz*(1-Al^2/2);%+Al^4/24);
   A2 = K2*Ma^2*Cm;
   C1 = K3*Ma^2*Cz;

   fg = [A1;A2;C1];
   fg = minlfr(fg,10000*eps);

% Now, the state-space matrices of the linearized models are computed
% by deriving the above expressions:

    ABCD = [diff(fg,'Al') diff(fg,'q') diff(fg,'dp')];

% Then, the input/output corresponding form is computed using
% "abcd2lfr" and the result is reduced using "minlfr".

   sys = abcd2lfr(ABCD,2);
   sys = minlfr(sys,10000*eps);

% The actuator is added in series at the system input:

   sys_nn = sys*ss(tf([omegaa^2],[1 2*kia*omegaa omegaa^2]));

% Finally, the the parameter variations are normalized.

   sys = normalizelfr(sys_nn);
   sys_6_2 = sys;

% The result of ex_6_1.m can be compared: distlfr(sys_6_1,sys_6_2)
