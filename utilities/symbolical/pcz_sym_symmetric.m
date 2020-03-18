function P = pcz_sym_symmetric(name,n,startindex)
%% pcz_sym_symmetric
%  
%  File: pcz_sym_symmetric.m
%  Directory: 2_demonstrations/lib/matlab
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2018. October 23.
%

%%

if nargin < 3
    startindex = 1;
end


nr = n*(n+1)/2;

p = pcz_sym_indexed(name, nr, startindex);

U = triu(ones(n));

P = sym(U);

P(U == 1) = p;

P = P + triu(P,1).';

end

function self_check
%%
% TODO: erre nem mukodik

P = pcz_sym_symmetric('p',4,0)

end