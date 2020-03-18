% LFR2ROB      - LFR-objects to 'USS' or 'UMAT' objects
%-----------------------------------------------------------------
% PURPOSE
% Generates a Robust Control Toolbox 3.0 Object from a LFR-object.
% Don't  support  inversion of objects with nominal value equal to
% zero.
%
% SYNOPSIS
% sys_out = lfr2rob(sys_in[,Ts]);
%
% INPUT ARGUMENTS
% sys_in  LFR-object.
% Ts      Sample time for discrete-time systems
%
% OUTPUT ARGUMENTS
% sys_out USS or UMAT object
%
% See also lfr2mu, lfr2mustab, lfr2mubnd, lfr2muss, lfr
% Robust Control Toolbox Version 3 required
%#----------------------------------------------------------------
% % EXAMPLE
%    lfrs x y z [1 2 3] [3 4 5]
%    M1 = [1+x^2*z (1+y)/(1+z^4) 1/(1-z);x^3 1/(2+y)^3 x^3/(2+y)^2];
%    size(M1)
%    M2 = lfr2rob(M1)
%    distlfr(M1,M2)
%
% % Check results
%    uplft(M1,'x',10,'y',20,'z',30)
%    usubs(M2,'x',10,'y',20,'z',30)
