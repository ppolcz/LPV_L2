function [ret] = pcz_plotpoint(V,varargin)
%% pcz_plotpoint
%  
%  File: pcz_plotpoint.m
%  Directory: utilities/plotting_tools
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2019. June 05. (2019a)
%

%%


[dim,N] = size(V);

assert(dim == 3 || dim == 2);

Z = zeros(1,N);

g = max(abs(V(:)));

if dim == 3
    plot3(V(1,:),V(2,:),V(3,:),'.',varargin{:});
else
    plot(V(1,:),V(2,:),'.',varargin{:});
end    

axis equal

if dim == 3
    axis([-g g -g g -g g])
else
    axis([-g g -g g])
end

end