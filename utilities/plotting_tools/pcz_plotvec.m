function [ret] = pcz_plotvec(V,varargin)
%% Script pcz_plotvec
%  
%  file:   pcz_plotvec.m
%  author: Peter Polcz <ppolcz@gmail.com> 
%  
%  Created on 2017. September 27.
%
%%

[dim,N] = size(V);

assert(dim == 3 || dim == 2);

Z = zeros(1,N);

g = max(abs(V(:)));

if dim == 3
    quiver3(Z,Z,Z,V(1,:),V(2,:),V(3,:),varargin{:});
else
    quiver(Z,Z,V(1,:),V(2,:),varargin{:});
end    

axis equal

if dim == 3
    axis([-g g -g g -g g])
else
    axis([-g g -g g])
end

end