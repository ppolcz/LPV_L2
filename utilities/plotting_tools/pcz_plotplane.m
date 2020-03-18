function [ret] = pcz_plotplane(V,varargin)
%% pcz_plotplane
%  
%  File: pcz_plotplane.m
%  Directory: 2_demonstrations/lib/matlab
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2018. September 17.
%

%%

[dim,N] = size(V);

assert(dim == 3 && N == 2);

V = orth(V);

g = 2*max(abs(V(:)));

[u,v] = meshgrid(linspace(-1,1,31));

surf(u*V(1,1)+v*V(1,2),u*V(2,1)+v*V(2,2),u*V(3,1)+v*V(3,2), varargin{:});

axis vis3d
axis([-g g -g g -g g])

end