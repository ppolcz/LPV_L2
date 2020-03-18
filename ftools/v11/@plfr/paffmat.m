function N = paffmat(pLFR,xi)
%% plfr2paffine
%  
%  File: plfr2paffine.m
%  Directory: 1_PhD_projects/24_passivity_ujsag/algebrai_probalkozasok_atalakitasok
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2019. November 06. (2019a)
%

%%

% Check whether, pLFR is indeed affine
if norm(pLFR.D) > 0
    error('Matrix D is not null, norm(D) = %d. LFR is not affine.', ...
        norm(pLFR.D))
end

% channels of the PAffineMatrix
if isa(xi,'lfr')
    xi_sym = sym(plfr(xi));
else
    xi_sym = sym(xi);
end
xi_sym = xi_sym(:);

nem_szerepel = setdiff(sym(pLFR.names),xi_sym);
if ~isempty(nem_szerepel)
    error('Variables `%s` must be included in xi!', ...
        strjoin(cellfun(@(x) {char(x)}, num2cell(nem_szerepel)),', '))
end

% Nr. of requested channels in the newly created PAffineMatrix.
s = numel(xi_sym);

% Size of each LFR blocks
r = pLFR.desc(1,:);
assert(size(pLFR.D,1) == sum(r), ...
    'Size not OK, not all blocks are taken into consideration')

% Helper matrices to produce kronecker product.
% e.g. [a,b] o [0 1 0 0] = [0 a 0 0 , 0 b 0 0]
E_cell = num2cell(eye(s),2);

% Find the indices of block names in xi (Eind).
tmp = sym(pLFR.names) ./ xi_sym;
tmp(tmp ~= 1) = 0;
[Eind,~] = find(tmp);

% Split up matrix C vertically corresponding to the block sizes
C_cell = cell(numel(r),1);
[C_cell{:}] = pcz_split_matrix(pLFR.C,r,Inf);

% pcz_symzero(vertcat(C_cell{:}) - pLFR.C, 'split of C')

% In each LFR block, the corresponding part of C should be multiplied as
% follows: Theta{i} = C{i} o [ 0 ... 1 ... 0 ], where 1 appears at index i,
% which denotes the index of desc.names{i} in xi.
Theta_cell = cellfun(@(C,i) {kron(C,E_cell{i})}, C_cell, num2cell(Eind));

Theta = vertcat(Theta_cell{:});
Im = eye(pLFR.nu);

N = PAffineMatrix(Theta,xi_sym);

% pcz_symzero_report(N.Sym - pLFR.Delta*pLFR.C, 'Theta * (Im o xi) = Delta * C')

N = pLFR.B*N + pLFR.A;

% pcz_symzero_report(N.Sym - sym(pLFR), 'pLFR = N');

end