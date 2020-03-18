% Illustration of "uplft" (feedback via the Delta-matrix)

  G = rlfr(4,2,3,5,6,'d_');

% Computation of the reduced LFR-object corresponding to
% s = 10j and delta_2 = 20:

   Gred1 = uplft(G,{'1/s','d_2'},[1/(10*j),20]);

% Compare

   size(G)
   size(Gred)

% Alternalively, using the function eval

   Int = 1/(10*j);
   d_2 = 20;
   Gred2 = eval(G);
