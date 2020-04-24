function [ret] = plus(X,Y)
%% plus
%  
%  File: plus.m
%  Directory: 1_PhD_projects/00_my_toolboxes/FinslerTools/v11/@PGenAffineMatrix
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2019. February 14.
%

%%

assert(all(size(X) == size(Y)), 'Size(X) = (%d,%d) != (%d,%d) = Size(Y)',...
    size(X), size(Y))

[X,Y] = adapt_channels(X,Y);

ret = PGenAffineMatrix(X.Theta + Y.Theta, X.channels, X.subsvars);
        
end


function test1
%%
syms p q r real
vars = [1;p;q;r];

X = PGenAffineMatrix(round(2*randn(4,12)),vars);
Y = PGenAffineMatrix(round(2*randn(4,12)),vars);

X.Sym = 1
Y.Sym = 1

Z = X + Y;

Z.Sym = 1

end

function test2
%%

tmp = PGenAffineMatrix(ones(3,6),p)
tmp + eye(3)
tmp = tmp + eye(3)
tmp.Sym = 1

end