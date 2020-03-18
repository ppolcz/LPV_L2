% NS_RAD       - Non-sigularity radius of an LFR-object
%-----------------------------------------------------------------
% PURPOSE
% Computation the non-sigularity radius of an lfr-object.
%
% SYNOPSIS
% [rad_min,rad_max,pertnames,pert] = ns_rad(A[,option]);
% ns_rad(A[,option]);
%
% DESCRIPTION
% The  bounds  (rad_min and rad_max) of the  radius  are  computed
% for  square  matrices.  For  rectangular matrices, there are two
% options for squaring:
%
% Option 0:  the  matrix A is squared down considering the nominal
% value A.d:  A -> [A;Q]  (where  Q is an orthogonal complement of
% A.d).  This  method induces artificial non-singularities, so, we
% obtain  a  bound  of the non-singularity radius the conservatism
% of which cannot be estimated.
%
% Option 1:   A  is  replaced by A*A'. In this case complex scalar
% uncertainties  are not supported and the size of the Delta block
% is  twice  the  original one (longer computation for mu). If the
% lower  bound  converges, we have an estimate of conservatism. In
% some  cases,  with option 0,  despite  the  possible  artificial
% singularities,  better results are obtained because the mu-upper
% bound has better convergence properties.
%
% INPUT ARGUMENTS
% A        lfr-object.
% option   0 or 1 (see above), option = 1 by default.
%
% OUTPUT ARGUMENTS
% rad_min  rad_min < non-sigularity radius < rad_max
% rad_max  We just have a guarantee that the considered LFR-object
%          remains  non-singular  for   norm(Delta) < rad_min. The
%          values  of  rad_min  and  rad_max permit us to evaluate
%          the conservatism.
% pert     If the lower bound of mu has converged (rad_max ~= inf)
%          pert  corresponds  to a possible the worst case. (Check
%          svd(A_wc.d) where A_wc =uplft(A,pertnames,pert).)
%
% See also wp_rad, lfr2mu, lfr2mustab, lfr2mubnd
%#----------------------------------------------------------------
% % EXAMPLE
%    lfrs a b
%
% % A1 is rank deficient for a = b = 0.4637
%    A1 = [1/(2+a) 1+a*b 2;-4/3+a+b -1-a*b -2];
%    ns_rad(A1);
%    ns_rad(A1,0); % here we are lucky to get a better result
%
% % Same with output arguments
%    [rad_min,rad_max,pertnames,pert] = ns_rad(A1,0);
%    X = uplft(A1,pertnames,pert);
%    % svd(X.d) rank shouldn't be maximal
