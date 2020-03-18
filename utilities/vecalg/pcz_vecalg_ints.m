function [Ints_X_Y] = pcz_vecalg_ints(X, Y, tol)
%% Script pcz_vecalg_ints
%  
%  File: pcz_vecalg_ints.m
%  Directory: projects/3_outsel/2018_01_10_LPV_inversion
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2018. January 13.
%
%%

if nargin < 3
    tol = 1e-10;
end

%%

Ker_X_UNION_Ker_Y = [pcz_vecalg_null(X',tol) , pcz_vecalg_null(Y',tol)];

Im__Ker_X_UNION_Ker_Y = pcz_vecalg_orth(Ker_X_UNION_Ker_Y,tol);

Ints_X_Y = pcz_vecalg_null(Im__Ker_X_UNION_Ker_Y',tol);

assert(pcz_vecalg_subset(Ints_X_Y,X,tol), ...
    '(X ∩ Y) ⊈ X (alongside the given rank tolerance: tol = %g.', tol);

assert(pcz_vecalg_subset(Ints_X_Y,Y,tol), ...
    '(X ∩ Y) ⊈ Y (alongside the given rank tolerance: tol = %g.', tol);

end