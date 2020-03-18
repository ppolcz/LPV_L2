function [indices, not_indices, N] = pcz_select_symvar(collection, selection)
%% pcz_select_symvar
%  
%  File: pcz_select_symvar.m
%  Directory: 2_demonstrations/lib/matlab
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2018. June 09.
%

%%

% collection
% selection

% collection = transpose(collection(:));
% selection = selection(selection(:));

s = numel(collection);
q = numel(selection);

if s < q
    [s,q] = deal(q,s);
    [collection,selection] = deal(selection,collection);
end

sel_fh = matlabFunction(selection, 'vars', {collection});

N = zeros(q,s);

% [1] generate the coefficients matrix
E = [1 ; zeros(s-1,1) ];
for i = 1:s
    N(:,i) = sel_fh(E);
    E = circshift(E,1);
end

indices = find(sum(N,1));
not_indices = setdiff(1:s,indices);

end