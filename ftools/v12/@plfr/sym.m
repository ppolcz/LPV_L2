function [pLFR_sym, PI_sym] = sym(pLFR)
%%
% Modified 2019.11.20. (november 20, szerda), 14:46

[A,B,C,D] = data(pLFR);  
[pLFR_sym, PI_sym] = sym_helper__(pLFR,A,B,C,D);

end
