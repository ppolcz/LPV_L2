function [Theta,z0,q] = P_Pi_canonical_decomp(z,w)
%% Script pcz_Pi_canonical_decomp
%  
%  file:   pcz_Pi_canonical_decomp.m
%  author: Peter Polcz <ppolcz@gmail.com> 
%  
%  Created on 2017.07.09. Sunday, 23:14:53
%
%%

TMP_wtNPjzxHKNoJIigzXrEl = pcz_dispFunctionName;

%%


if nargin == 1
    w = symvar(z);
end

m = numel(z);
r = sym('NONAMETMPVAR',[1 m]);

[zp,q] = numden(simplify(r*z));

[qc,~] = coeffs(q,w);
q = q / qc(1);

[c,t] = coeffs(expand(zp/qc(1)),w);
[A,~] = equationsToMatrix(c,r);

Theta = double(A');
z0 = t';

%%
pcz_dispFunctionEnd(TMP_wtNPjzxHKNoJIigzXrEl);
clear TMP_*

end

