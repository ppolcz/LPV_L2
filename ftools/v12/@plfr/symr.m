function [pLFR_sym, PI_sym, fancy, s] = symr(pLFR,varargin)
%%
%  File: symr.m
%  Directory: 7_ftools/ftools/v12/@plfr
%  Author: Peter Polcz (ppolcz@gmail.com) 
% 
%  Created on 2019.11.20. (november 20, szerda), 14:46
%  Major review on 2020. April 10. (2019b)
% 

opts.tol = 1e-10;
opts.maxden = 1000;
opts.maxit = 1;
opts.structout = true;
opts.sparsity = 15; % percent
opts = parsepropval(opts, varargin{:});

pLFR_sym = [];
PI_sym = [];
fancy = []; 
s = [];

nonzero_elemets_nr = numel(find(abs(pLFR.M) > opts.tol));
nonzero_elemets_perc = nonzero_elemets_nr / numel(pLFR.M) * 100;

if nonzero_elemets_nr > opts.sparsity
    return
end

r = @(M) pcz_find_recdec(M,'maxden',opts.maxden,'maxit',opts.maxit,'tol',opts.tol, 'structout', opts.structout);

[M,s] = r(pLFR.M); 

fancy = s.good;

if nargin > 1 && ~fancy
    return
end

if fancy
    [A,B,C,D] = pcz_split_matrix(M,[Inf pLFR.m1],[Inf pLFR.m1]);
else
    [A,B,C,D] = data(pLFR);
end

[pLFR_sym, PI_sym] = sym_helper__(pLFR,A,B,C,D);

end
