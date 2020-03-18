% LFR2MUSTAB   - LFR-objects adapted to the function mustab
%-----------------------------------------------------------------
% PURPOSE
% Generates  the  input arguments of the function 'mustab' for mu-
% analysis of a given LFR-object.
%
% SYNOPSIS
% [P,delta] = lfr2lmip(sys[,K]);
%
% REMARK
% With version 3 of the Robust Control Toolbox you may try
% robuststab(lfr2rob(sys))
%
% INPUT ARGUMENTS
% sys     Dynamic LFR-object.
% K       Feedback:  matrix, ss-object, lfr-object (scheduled fb).
%
% OUTPUT ARGUMENTS
% P       LTI-system (compatible with pck, unpck).
% delta   Uncertainty  block  description  compatible with the LMI
%         Control toolbox.
%
% See also lfr2mu, lfr2mubnd, lfr2mussv
%#----------------------------------------------------------------
% % EXAMPLE
% % Uncertain system
%    lfrs a b c
%    A = [-1+5*a -10*(1+b);10*(1+c^2) -1-a*b*c];
%    B = [1+b;1*c]; C = [1 1];
%    sys = abcd2lfr([A B;C 0],2);
%
% % Computation of the arguments of mubnd
%    [P,delta] = lfr2mustab(sys);
%    [margin,frequ0] = mustab(P,delta);
%
% % Result (For comparison, see help lfr2mubnd)
%    mu_mx = 1/margin
