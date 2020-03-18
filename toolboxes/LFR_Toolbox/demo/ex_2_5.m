% Illustration of lower LFT (standard feedback but which removes
% inputs/outputs)

   G = rlfr(4,2,3,5,6,'d_');
   K = rand(3,2);

% Feedback

   Gfb = zeros(0,2)*feedback(G,K,1)*zeros(3,0);

% Results

   size(Gfb)

% The new system has no longer inputs and outputs
