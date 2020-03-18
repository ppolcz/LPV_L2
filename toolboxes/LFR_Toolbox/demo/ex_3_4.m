% Illustration of the function gmorton (generalized Morton's
% realization) and of the function abcd2lfr.

% Definition of the coefficients of an expansion

   sys0 = rss(4,2,3);
   sys1 = rss(4,2,3);
   sys2 = rss(4,2,3);

% Realization of the system with [A B;C D] of the form
% [A B;C D] = [A0 B0;C0 D0] + d1*[A1 B1;C1 D1] + d2*[A2 B2;C2 D2]

   lfrs d1 d2
   sys = gmorton({sys0,sys1,sys2},[1 d1 d2]);

% For checking the result

   abcd = [sys0.a sys0.b;sys0.c sys0.d] + ...
   d1 * [sys1.a sys1.b;sys1.c sys1.d] + ...
   d2 * [sys2.a sys2.b;sys2.c sys2.d];

   newsys = abcd2lfr(abcd,4);

   distlfr(sys,newsys)

% This function is a generalization of the Morton's method because
% we can use complex expansion form:

   lfrs d1 d2
   renewsys = gmorton({sys0,sys1,sys2},[1+d1^2 d1*d2 1/(1+d2)]);


