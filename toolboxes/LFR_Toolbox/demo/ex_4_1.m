% This example illustrates order reduction after realization (here
% object-oriented realization) using the 1-D approach. The fact that
% parameter commutativity is ignored is also illustrated.

% Realization:

   lfrs d1 d2
   S = [d1*d2;d1*d2];

% This object is equivalent to [1;1]*d1*d2 i.e. of minimum order 2.

% The function for 1-D reduction is "minlfr1". By default it reduces
% the order, first w.r.t. 1/s (empty block here), and then w.r.t.
% the first, the second (and so on) uncertain parameter.

   Smin = minlfr1(S);

% resulting in

   size(Smin)

% Reduction with respect to the first uncertain parameter is not
% found because d1 cannot be factorized on the right (commutativity
% is ignored by "minlfr1").

% In order to apply 1-D reductions more than one time for each
% parameter we must use an optional parameter ([2 3 2] below):

   Smin2 = minlfr1(S,[2 3 2]);

% In the second argument, ordering starts from 1/s so, "2" stands
% for the first uncertain parameter, "3" for the second one and so on.
% So, 1-D reduction is performed respectively relative to d1, d2 and d1.
% The result is:

   size(Smin2)

% As expected, now, both factorizations are "found"
