function pLFR = generate_dimensions(pLFR)
%%

% LFR dimensions (output, input)
[pLFR.ny,pLFR.nu] = size(pLFR.A);

% LFR dimensions (Delta block)
pLFR.m1 = size(pLFR.Delta,1);

end
