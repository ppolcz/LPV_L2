% Loop invoked by ex_6_5.m

   kk = 1;

   warning off
   for ii = 0:5;
   for jj = 0:5;

       % Trim parameters
       Al = ii*(0.349/5); % from 0 to 0.349
       Ma = 2 + jj*(2/5); % from 2 to 4

       % Computation of equilibrium points
       X0 = [0;0;Al;0];
       [X,U,Y,dX] = trim('missile53',X0,0,0,3);

       % Computation of linearized models
       [A,B,C,D] = linmod('missile53',X,U);
       B = B/B(1,1);
       C = C*B(1,1);

       % Storage of values for state-space model interpolation
       abcd_data{kk} = [A B;C D];
       al_ma_data{kk} = [Al,Ma];

       % Storage of values for equilibrium surface interpolation
       q_dp_data{kk} = [X(4) U];  % values of q and dp

       kk = kk + 1;

   end;
   end;
   warning on
