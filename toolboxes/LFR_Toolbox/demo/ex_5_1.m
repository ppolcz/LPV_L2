% This example illustrates modelling of a full dynamic uncertainty
% block

% The function lfr offers an option for defining full blocks with
% frequency domain bounds. However, this feature is not already
% interfaced with MATLAB functions, therefore, it is suggested
% to build normalized full blocks and to add manually the weighting
% function W(s).
%
%                 +---[W]---[Delta]---+
%                 |                   |
%          ------>o-------------------+---[M]--->

% An lfr-system
   M = rlfr(5,2,3,4,4,4,'m');

% Weightings
   W = ss(tf([30 100],[1 100])^2)*eye(2,2);

% Full complex block
   Delta = lfr('Delta','ltifc',[2,2]);

% Now we put in series M and the identity with feedback Delta*W
   sys = feedback(lfr(eye(2,2)),Delta*W,1)*M;
   size(sys)

