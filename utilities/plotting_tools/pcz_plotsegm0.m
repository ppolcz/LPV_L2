function [ret] = pcz_plotsegm0(V,varargin)
%% pcz_plotsegm(V,varargin)
%  
%  File: pcz_plotsegm(V,varargin).m
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
    P = plot3([0 V(1,1)],[0 V(2,1)],[0 V(3,1)],'.-',varargin{:});
else
    P = plot([0 V(1,1)],[0 V(2,1)],'.-',varargin{:});
end
hold on

for i = 2:size(V,2)
    if dim == 3
        plot3([0 V(1,i)],[0 V(2,i)],[0 V(3,i)],'.-','Color',P.Color,varargin{:});
    else
        plot([0 V(1,i)],[0 V(2,i)],'.-','Color',P.Color,varargin{:});
    end    
end

axis equal

if dim == 3
    axis([-g g -g g -g g])
else
    axis([-g g -g g])
end

end