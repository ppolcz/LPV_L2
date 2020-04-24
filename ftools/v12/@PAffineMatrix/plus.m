function [ret] = plus(X,Y)
%% plus
%  
%  File: plus.m
%  Directory: 1_PhD_projects/00_my_toolboxes/FinslerTools/v11/@PAffineMatrix
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2019. February 14.
%

%%

assert(all(size(X) == size(Y)), 'Size(X) = (%d,%d) != (%d,%d) = Size(Y)',...
    size(X), size(Y))

[X,Y] = adapt_channels(X,Y);

ret = PAffineMatrix(X.Theta + Y.Theta, X.channels, X.subsvars);
        
end

