% LFR2MUBND    - LFR-objects adapted to the function mubnd
%-----------------------------------------------------------------
% PURPOSE
% Generates  the  input  arguments of the function 'mubdn' for mu-
% analysis of a given LFR-object.
%
% SYNOPSIS
% [M,delta] = lfr2mubnd(sys,frequ[,K]);
%
% INPUT ARGUMENTS
% sys     LFR-object
% frequ   Real number = frequency at which M must be computed.
%         frequ can be a vector, in this case M is a cell.
% K       Feedback:  matrix, ss-object, lfr-object (scheduled fb).
%
% OUTPUT ARGUMENTS
% M       Constant  matrix (the M-matrix in an M-Delta form). Cell
%         containing such matrices if frequ is a vector.
% delta   Uncertainty  block  description  compatible with the LMI
%         Control toolbox.
%
% See also lfr2mustab, lfr2mu, lfr2mussv
%#----------------------------------------------------------------
% % EXAMPLE
% % Uncertain system
%    lfrs a b c
%    A = [-1+5*a -10*(1+b);10*(1+c^2) -1-a*b*c];
%    B = [1+b;1*c]; C = [1 1];
%    sys = abcd2lfr([A B;C 0],2);
%
% % Computation of the arguments of mubnd
%    frequ = logspace(0,2,20);
%    [M,delta] = lfr2mubnd(sys,frequ);
%
% % Mu-analysis
%    bnds = [];
%    for ii = 1:length(frequ);
%       bnds(ii) = mubnd(M{ii},delta,1e-6,'off');
%    end;
%
% % Results (For comparison, see help lfr2mu)
%    semilogx(frequ,bnds,'g-'); hold on;
%    mu_max = max(bnds); % See help lfr2mustab
