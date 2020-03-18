% DATA2SYM     - Interpolation by mean squares
%-----------------------------------------------------------------
% PURPOSE
% Mean squares interpolation of matrices (e.g., matrices [A B;C D]
% of  a linear system), on a basis of symbolic or lfr-objects. The
% result is a symbolic object.
%
% SYNOPSIS
% [Symform,Kcell,Kfix,Err_min,Err_max,Kerror] =
%        data2sym(Kdata,pardata,symex,order[,tk[,Kcstr,parcstr]]);
%
% DESCRIPTION
% See  data2lfr. This function is very similar to data2lfr but the
% result  is  in  a  symbolic form for permitting the user to call
% symtreed for reduced order lfr-object realization.
%
% INPUT ARGUMENTS
% Kdata    Cell  containing the matrices to be interpolated.
% pardata  Cell  containing the vectors of the parameter values at
%          interpolation points.
% symex    Vector  of  (1-by-1)  symbolic or lfr-objects involving
%          the objects listed in 'order' (interpolation formula).
% order    Cell of 1-by-1 symbolic or lfr-objects for ordering the
%          numerical data contained in pardata{ii}:
%          pardata{ii}(jj) = numerical  value of object order{jj}.
% tk       Talkative indice: tk = 0,1 or 2.  No message if tk = 0.
% Kcstr    Optional.
% parcstr  Optional.  This  pair of arguments is similar to a pair
%          of  entries  of Kdata and pardata: defines a constraint
%          to be exactly met (no least squares).
%
% OUTPUT ARGUMENTS
% Symform  lfr-object interpolating data in 'Kdata' and 'pardata'.
% Kcell    Coefficients of the interpolation formula ({K0,K1...}).
% Kfix     Entries that do not vary in K0, K1, K2....
% Err_min  Minimum absolute error entry per entry.
% Err_max  Maximum absolute error entry per entry.
%          Err_min  and  Err_max can be used for modelling uncert-
%          ainties due to interpolation (see bnds2lfr).
% Kerror   Cell  containing  interpolation relative errors at each
%          point.
%
% See also data2lfr, symtreed, bnds2lfr, abcd2lfr
%#----------------------------------------------------------------
% % EXAMPLE
% % Data
%     K=[]; K{1}=rand(2,3); K{2}=rand(2,3); K{3}=rand(2,3);
%     p=[]; p{1}=rand(1,2); p{2}=rand(1,2); p{3}=rand(1,2);
%     Kcstr = rand(2,3);
%
% % Interpolation formula K0 + a K1 + b K2 + a b K3
% % 1- Interpolation without constraints
%     syms a b
%     Symform1  = data2sym(K,p,[1 a b],{b a},2);
%
% % 2- Interpolation with constraints
%     lfrs a b
%     Symform2 = data2sym(K,p,[1 a b a*b],{a b},2,Kcstr,[4 7]);
%
% % Verification of the constraint
%     Kcstr2 = uplft(Lfrform2,{a b},[4 7]);
%     Kcstr2.d - Kcstr
