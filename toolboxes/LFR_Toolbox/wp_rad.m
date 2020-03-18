% WP_RAD       - Well-posedness radius of an LFR-object
%-----------------------------------------------------------------
% PURPOSE
% Computes the well-posedness radius of an LFR-object
%
% SYNOPSIS
% [rad_min,rad_max,pertnames,pert] = wp_rad(A);
% wp_rad(A);
%
% INPUT ARGUMENTS
% A        lfr-object
%
% OUTPUT ARGUMENTS
% rad_min  rad_min < well-posedness radius < rad_max
% rad_max  We just have a guarantee that the considered LFR-object
%          remains  well-posed  for  norm(Delta)  <  rad_min.  The
%          difference  between  the  values of rad_min and rad_max
%          permit us to evaluate the conservatism of the test.
% pert     If the lower bound of mu has converged (rad_max ~= inf)
%          pert corresponds to a worst case.
%
% See also ns_rad, lfr2mu, lfr2mustab, lfr2mubnd
%#----------------------------------------------------------------
% % EXAMPLE
%    lfrs a b
% % A is not well posed for a = -2
%    A = [1/(2+a) 1+a*b 2;-4/(4.1+a+b) -1-a*b -2];
%    wp_rad(A)
