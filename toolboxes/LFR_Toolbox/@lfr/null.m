% NULL         - Kernel of an LFR-object
%-----------------------------------------------------------------
% PURPOSE
% Computes  the  kernel  of  a lfr-object using a simple and fast
% methode.  The  counterpart of fastness is that the radius of the
% ball in which the result is well-posed is not optimized.
%
% SYNOPSIS
% X = null(A[,tol])
%
% DESCRIPTION
% The  lfr-object  A  is  squared down using the matrix denoted Q
% below. The kernel is  computed  as follows:
%              X = inv([A;Q])*[0;I].
% The  main  problem of this approach is that the row span of Q is
% likely  to  intersect  the row span of A for some combination of
% parameters.  The matrix [A;Q] cannot be inverted at these points
% reducing the well-poseness radius of the computed kernel.
%
% INPUT ARGUMENTS
% A     lfr-object such that A.d has maximum row rank.
% tol   Tolerance for reduction (by minlfr).
%
% See also wp_rad2, ns_rad2
%#----------------------------------------------------------------
% % EXAMPLE
%    lfrs a b
%    A = [1+a*b^2 1+b^3 2;0 b 3+a*(1+b^2)];
%    X = null(A);
%    size(X)
%
%    distlfr(A*X,zeros(2,1))
