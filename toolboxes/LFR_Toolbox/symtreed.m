% SYMTREED     - Elementary structured tree decomposition
%-----------------------------------------------------------------
% PURPOSE
% Computes  the lfr-object corresponding to a polynomial matrix in
% symbolic form using the tree decomposition approach for low oder
% modelling.
%
% SYNOPSIS
% lfrout = symtreed(symbex[,reorder[,isverbose]])
%
% DESCRIPTION
% The  technique  used here is described in J.C. Cockburn and B.G.
% Morton,  Automatica  Vol. 33, No 7, pp 1263-1271, 1997.  The low
% order  model  obtained  by this technique can be further reduced
% by applying minlfr or minlfr1.
%
% NOTE
% - The denominators of the first input argument must be constant.
% See  the  manual  for  explanations concerning the fact that all
% rational  objects  can be transformed to polynomial objects (see
% also help lf2lfr or help rf2lfr).
% - The  main  limitation  is the following: If the input argument
% symbex  is  given  as a complex factorized formula, it is likely
% that  sym2lfr (which doesn't remove factorizations) will lead to
% much  better  results  than symtreed (which expands numerators).
% For  example  M = (1+a^2+b*a)^5-(1+a+b)^3 leads to order 26 with
% sym2lfr,  and  to  order  39  with  symtreed (note  that if M is
% replaced  by  expand(M)  before  sym2lfr  is  invoked, the order
% becomes 156, it means that symtreed has reduced form 156 to 39).
%
% INPUT ARGUMENTS
% symbex      Symbolic POLYNOMIAL matrix.
% reorder     [] or cell of strings of length less or equal to the
%             number  of  parameters of symbex. These strings must
%             be parameter names. This argument is used for reord-
%             ering  the  parameters during the tree decomposition
%             process. This feature is important in some cases for
%             obtaining lower order lfr-objects.
% isverbose   0 (default) or 1.  If equal to zero, no comments are
%             displayed.
%
% OUTPUT ARGUMENT
% lfrout      Returned lfr-object.
%
% See also sym2lfr
%#----------------------------------------------------------------
% % EXAMPLE
%    syms Int x y z
%    symbex = [Int*x^2*z^3+3*x*y*z,x^2*y^2*z^3,1;0,Int^2*x^3*y*z,x*y^2*z-x^2*y*z]
%
% % Direct realization
%    sys0 = sym2lfr(symbex,-1);
%
% % Realization using an elementary structured tree decomposition
%    sys1 = symtreed(symbex);
%
% % Similar with reordering of parameters
%    sys2 = symtreed(symbex,{'z','y','Int','x'});
%
% % Results:
%    size(sys0)
%    size(sys1)
%    size(sys2)
