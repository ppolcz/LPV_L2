function [lfr] = P_lfrdata_v6(G_lfr)
%% Script P_lfrdata_v6
%  
%  file:   P_lfrdata_v6.m
%  author: Peter Polcz <ppolcz@gmail.com> 
%  
%  Created on 2017.08.01. Tuesday, 00:15:12
%
%%

TMP_EIHwhNmHSeYbRRcCXqPF = pcz_dispFunctionName;

%%


[D,C,B,A,blk] = lfrdata(G_lfr);
names = blk.names;
desc = blk.desc(1:2,:)';

% Delta
s = numel(names);
c = cell(s, 1);
for k = 1:s
    c{k} = sym(names{k}) * eye(desc(k,:));
end
Delta = blkdiag(c{:});

[m,~] = size(C);
I = eye(m);

% G(x)x + F(x)π = Cb * πb = 0
G = -Delta*C;
F = I-Delta*D;

lfr = v2struct(A,B,C,D,I,blk,Delta,G,F);

%%
% End of the script.
pcz_dispFunctionEnd(TMP_EIHwhNmHSeYbRRcCXqPF);
clear TMP_*

end


