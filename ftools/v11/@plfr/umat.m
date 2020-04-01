function [G_umat] = umat(G)
%% umat
%  
%  File: umat.m
%  Directory: 7_ftools/ftools/v11/@plfr
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2020. March 30. (2019b)
%

%%

deltai_fh_cell = cellfun(@(s) { matlabFunction(s,'vars',G.symvars) }, num2cell(diag(G.Delta)));

n = cell(G.varnames).';
b = num2cell(G.bounds,2);
% r = num2cell(G.rbounds,2);

ureal_cell = cellfun(@(n,b) { ureal(n,sum(b)/2,'Range',b) }, n, b);

deltai_ureal_cell = cellfun(@(d) { d(ureal_cell{:}) }, deltai_fh_cell);

Delta_umat = blkdiag(deltai_ureal_cell{:});

G_umat = G.A + G.B * ( (G.I - Delta_umat * G.D) \ Delta_umat ) * G.C;

end