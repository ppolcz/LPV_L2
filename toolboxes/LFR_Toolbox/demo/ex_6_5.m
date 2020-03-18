% Example 6.5: trim / linearization on a gridding. Interpolation
% for equilibrium surface and for state-space matrices.

   missiledata

   % Trim, linearization and storage of results at gridding points
   ex_6_5_loop;

% Results of ex_6_5_loop
% al_ma_data = values of alpha and Mach defining the gridding
% q_dp_data  = trimmed values of q and dq at gridding points

% ================================================================
% Interpolation for finding the equilibrium surface
% It means finding q and dp as a polynomial expansion of Al and Ma

   lfrs Al Ma [0.0 2.0] [0.349 4.0]

   lfrex = [1 Al Al*Ma Al^2];
   ordlfr = {Al Ma};
   q_dp = data2lfr(q_dp_data,al_ma_data,lfrex,ordlfr);

   % q at equilibrium as a function of alpha and Mach
   q = q_dp(1);

   % dp at equilibrium as a function of alpha and Mach
   dp = q_dp(2);

   % Verification of results (dp2 is the exact equation)
   dp2 = -(m3*Al^3 + m2*Al^2 + m1*(-7 +(8/3)*Ma)*Al)/m0;
   distlfr(dp,dp2)

   % Verification by plotting the surfaces dp and dp2
   figure
   plotlfr(dp, {'Al',0,0.349,10},{'Ma',2,4,10});
   hold on
   plotlfr(dp2,{'Al',0,0.349,10},{'Ma',2,4,10});

