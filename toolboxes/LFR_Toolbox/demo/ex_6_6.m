% Example 6.6: trim / linearization on a gridding. Interpolation
% for state-space matrices.

   missiledata

   % Trim, linearization and storage of results at gridding points
   ex_6_5_loop;

% Results of ex_6_5_loop
% al_ma_data = values of alpha and Mach defining the gridding
% abcd_data  = state-space matrices of the linearized models at
%              gridding points.

% ================================================================
% Interpolation for finding the continuum of linearized models
% It means finding [A B;C D] as a polynomial expansion of Al and Ma

   lfrs Al Ma [0.0 2.0] [0.349 4.0]

   % Interpolation
   lfrex = [1 Al Ma Al*Ma Al^2 Ma^2 Al^2*Ma Al*Ma^2 Al^3 Ma^3];
   ordlfr = {Al Ma};
   ABCD = data2lfr(abcd_data,al_ma_data,lfrex,ordlfr);

   % Order reduction
   ABCD = minlfr(ABCD,0.00000001);
   size(ABCD)

   % From (A,B,C,D) to input/output
   sys_nn = abcd2lfr(ABCD,4);

   % Normalization
   sys = normalizelfr(sys_nn);

   sys_6_6 = sys;
   size(sys_6_6)

