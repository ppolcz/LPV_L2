% Illustration of the function lfr/eval (see also ex_2_7.m)
% Advanced "star product"

% Random lfr-object with parameters d1 and d2

   M1 = rlfr(4,2,3,4,4,'d');
   size(M1)

% Definition of d1 and d2 in the workspace
% Equivalent to lfrs d1 d2

   d1 = lfr('d1','ltisr');
   d2 = lfr('d2','ltisr');

% Define d2 as a function of d1

   d2 = d1^2;

% Evaluation

   M2 = eval(M1);
   M2.blk.names

% d2 vanished

% Alternative evaluation replacing 1/s by a function of 1/z
% (Tustin's transformation)

   Delay = lfr(0,1,1,0,'d');
   Int =0.1 * (1+Delay) / (1-Delay);

% Evaluation

   M3  = eval(M2);
   M3.blk.names
