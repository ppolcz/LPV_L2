% DATA2LFR     - Interpolation by mean squares
%-----------------------------------------------------------------
% PURPOSE
% Mean squares interpolation of matrices (e.g., matrices [A B;C D]
% of a linear system), on a basis of lfr-objects. The result is an
% lfr-object.
%
% SYNOPSIS
% [Lfrform,Kcell,Kfix,Err_min,Err_max,Kerror] =
%       data2lfr(Kdata,pardata,lfrex,order[,tk[,Kcstr,parcstr]]);
%
% EXAMPLE
% We  look  for  the  matrices  K0, K1, K2 such that the following
% equality is satisfied (least squares sense):
% K(a,b)   = 1*K0 + a*K1 + 3*a^2*b*K2
% Kdata    = { K(a1,b1) , K(a2,b2),... }  contains a set of numer-
%            ical values of K(a,b).
% pardata  = { [a1,b1] , [a2,b2],... } contains a set of numerical
%            values of a and b.
% lfrex    = [1 a a^2*b]  contains  the  basis  of lfr-objects for
%            interpolation.
% order    = {a,b}  contains  a list of 1-by-1 lfr-objects ordered
%            like in 'pardata'.
% Results:  Kcell = {K0 , K1,...} and
% Lfrform = Kfix + lfrex(1) Kcell{1} + lfrex(2) Kcell{2} + ...
%
% This function does not attempt to interpolate the entries of the
% matrices  contained  in  Kdata  that  are  invariant  w.r.t. the
% interpolation points (-> Kfix).
%
% CAUTION
% If  some  numerical  values in 'pardata' are much larger than 1,
% you MUST scale them (e.g. replaces meters by kilometers)
%
% INPUT ARGUMENTS
% Kdata    Cell  containing the matrices to be interpolated.
% pardata  Cell  containing the vectors of the parameter values at
%          interpolation points.
% lfrex    Vector  of  (1-by-1) lfr-expressions involving the lfr-
%          objects  listed  in 'order' (interpolation formula).
% order    Cell  of 1-by-1 lfr-objects for  ordering the numerical
%          data contained in pardata{ii} ->
%          pardata{ii}(jj) = numerical value  of object order{jj}.
% tk       Talkative indice: tk = 0,1 or 2.  No message if tk = 0.
% Kcstr    Optional.
% parcstr  Optional.  This  pair of arguments is similar to a pair
%          of  entries  of Kdata and pardata: defines a constraint
%          to be exactly met (no least squares).
%
% OUTPUT ARGUMENTS
% Lfrform  lfr-object interpolating data in 'Kdata' and 'pardata'.
% Kcell    Coefficients of the interpolation formula ({K0,K1...}).
% Kfix     Entries that do not vary in K0, K1, K2....
% Err_min  Minimum absolute error entry per entry.
% Err_max  Maximum absolute error entry per entry.
%          Err_min  and  Err_max can be used for modelling uncert-
%          ainties due to interpolation (see bnds2lfr).
% Kerror   Cell  containing  interpolation relative errors at each
%          point.
%
% See also data2sym, bnds2lfr, abcd2lfr
%#----------------------------------------------------------------
% % EXAMPLE
% % Data
%     K=[]; K{1}=rand(2,3); K{2}=rand(2,3); K{3}=rand(2,3);
%     p=[]; p{1}=rand(1,2); p{2}=rand(1,2); p{3}=rand(1,2);
%     Kcstr = rand(2,3);
%
% % Interpolation formula K0 + a K1 + b K2 + a b K3
% % 1- Interpolation without constraints
%     lfrs a b
%     Lfrform1  = data2lfr(K,p,[1 a b],{b a},2);
%
% % 2- Interpolation with constraints
%     lfrs a b
%     Lfrform2 = data2lfr(K,p,[1 a b a*b],{a b},2,Kcstr,[4 7]);
%
% % Verification of the constraint
%     Kcstr2 = uplft(Lfrform2,{a b},[4 7]);
%     Kcstr2.d - Kcstr
