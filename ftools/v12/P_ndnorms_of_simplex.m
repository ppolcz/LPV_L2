function [v, fci, proj_x, proj_xp, in_simplex, Vk] = P_ndnorms_of_simplex(val,offset,p_lim,varargin)
%% P_ndnorms_of_simplex
%  
%  File: P_ndnorms_of_simplex.m
%  Directory: 7_ftools/ftools/v11
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2019. November 01. (2019a)
%

if nargin < 2
    offset = 0;
end

%%

nx = numel(val);
z = zeros(nx,1);

np = 0;
if nargin >= 3
    np = size(p_lim,1);
end

if isscalar(offset)
    offset = offset*(z+1);
elseif numel(offset) ~= nx
    error('dim(vals) ~= dim(offset)')
end

p = offset(:)' + z;
p = p + diag(val);

v = [
    offset(:)'
    p
    ];

fci = combnk(2:nx+1,nx-1);

fci = [
    ones(size(fci,1),1) fci
    2:nx+1
    ];
    
Nr = size(fci,1);
Vk = cell(Nr,1);
proj_x = cell(Nr,1);
proj_xp = cell(Nr,1);

for i = 1:Nr

    r0 = v(fci(i,1), :)';
    V = (  v(fci(i,2:nx),:) - repmat(r0', [nx-1, 1])  )';
        
    % normal vector (norm(n) == 1)
    nvec = null(V');

    % distance of the origin to the dim dimensional facet
    dst = dot(nvec, r0);

    % [3.3] Vk
    V = orth(V);
    Vk{i} = V;
    
    % [3.4.1] projection (only x needed to be substituted)
    % Renamings: 2019.11.21. (november 21, csütörtök), 13:43
    proj_x{i} = @(x) V*V'*x + dst*nvec*ones(1,size(x,2));
    
    % 2019.11.21. (november 21, csütörtök), 13:43
    V_xp = blkdiag(V,eye(np));
    nvec_xp = [nvec ; zeros(np,1)];

    % [3.4.2] projection in extended space ([x;p] needed to be substituted)
    % 2019.11.21. (november 21, csütörtök), 13:43
    proj_xp{i} = @(xp) V_xp*V_xp'*xp + dst*nvec_xp*ones(1,size(xp,2));
end


V = Vk{Nr};
nvec = null(V');
dst = dot(nvec', p(:,1)');

in_simplex = @(x) nvec'*x <= dst + 1e-10;



end