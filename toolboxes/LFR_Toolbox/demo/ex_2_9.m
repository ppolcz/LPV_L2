% Illustration of the "constant block" introduced for "inverting"
% non-invertible LFR-objects. Normalization and converse.

   S = rlfr(5,3,3,2,2,'d');
   S.d = zeros(3,3);

% "Inversion" of the non-invertible LFR-object

   invS = inv(S);
   size(invS)

% The second line displayed by the function size shows that the
% dummy parameter is repeated 3 times in a block named "constant
% block".

% Modification of the bound information. Here, the nominal value (4)
% is not centered in the range of variations ([2 8])

  set(invS,{'d1','d2'},{[2 8 4],[2 8 4]},{'minmax','minmax'});
  size(invS)

% The function 'set' changes the bound information stored in the
% LFR-object but bo not modify numerical values. For updating the
% nominal value, the function normalizelfr must be invoked.

  invS2 = normalizelfr(invS);

% Inversion is feasible now, therefore the dummy parameter has
% disappeared as shown below

   size(invS2)

% Note that it is possible to compute the actual values of parameters
% from normalized ones.

   par_nom = [-0.5 0.5];
   par_act = actualval(invS2,{'d1','d2'},par_nom)

% uplft(invS,{'d1','d2'},par_act) and uplft(invS2,{'d1','d2'},par_nom)
% are equal.

% The system invS2 can be unnormalized

   distlfr(invS,unnormalize(invS2))
