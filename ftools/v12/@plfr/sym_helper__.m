function [pLFR_sym, PI_sym] = sym_helper__(pLFR,A,B,C,D)
%%
% 2019.11.20. (november 20, szerda), 14:46

PI_sym = [
    eye(pLFR.nu)
    (pLFR.I-pLFR.Delta*D)\pLFR.Delta*C 
    ];

if ~isnumeric(PI_sym)
    PI_sym = simplify(PI_sym);
end

pLFR_sym = [A B]*PI_sym;

if ~isnumeric(pLFR_sym)
    pLFR_sym = simplify(pLFR_sym);
end

end