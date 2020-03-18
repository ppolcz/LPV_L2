% Illustration of the functions lfr/set and normalizelfr

% First a random LFR-object with parameters d1,d2,d3,d4 is defined

   S = rlfr(4,2,2,4,1,1,1,'d');
   size(S)

% Bounds are defined for the four parameters using the function
% set

   set(S,{'d1'},{[-2,6]},{'minmax'});
   set(S,{'d2'},{[-6,2]},{'minmax'});
   set(S,{'d3'},{[-2,2]},{'minmax'});
   set(S,{'d4'},{[-1,1]},{'minmax'});

   size(S)

% Normamalization

   Snorm = normalizelfr(S,{'d1','d2','d3','d4'});
   size(Snorm)

% For checking the results, both LFR-objects are evaluated on the
% same vertex

   uplft(S,{'d1','d2','d3','d4'},[6,-6,2,-1])
   uplft(Snorm,{'d1','d2','d3','d4'},[1,-1,1,-1])
