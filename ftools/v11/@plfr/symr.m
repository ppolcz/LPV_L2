function [pLFR_sym, PI_sym, fancy, s] = symr(pLFR,~)
%%
% Created: 2019.11.20. (november 20, szerda), 14:46

[A,B,C,D] = data(pLFR);

r = @(M) pcz_find_recdec(M,'maxden',1000,'maxit',1,'tol',1e-10, 'structout', true);

[M,s] = r(pLFR.M); 

fancy = s.good;

if nargin > 1 && ~fancy
    pLFR_sym = [];
    PI_sym = [];
    return
end

if fancy
    [A,B,C,D] = pcz_split_matrix(M,[Inf pLFR.m1],[Inf pLFR.m1]);
end

[pLFR_sym, PI_sym] = sym_helper__(pLFR,A,B,C,D);

end
