function [ret] = mtimes(X,Y)
%% mtimes
%  
%  File: mtimes.m
%  Directory: 1_PhD_projects/00_my_toolboxes/algo_P/@PGenAffineMatrix
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2018. October 21.
%

%%

XisAffinMat = isa(X,'PGenAffineMatrix');
YisAffinMat = isa(Y,'PGenAffineMatrix');


if XisAffinMat && ~YisAffinMat
    
    ret_Theta = X.Theta * kron(Y,X.Is);
    
    ret = PGenAffineMatrix(ret_Theta,X.channels,X.subsvars);
    
    if X.issym, ret.Sym = 1; end
    
    return

elseif ~XisAffinMat && YisAffinMat
    
    ret_Theta = X * Y.Theta;
    
    ret = PGenAffineMatrix(ret_Theta,Y.channels,Y.subsvars);
    
    if Y.issym, ret.Sym = 1; end
    
    return
end


ret = [];

end


function test1
%%

Theta = [
    1 0 0 0 0 1
    0 1 0 1 0 0
    ];

syms p q real

A = PGenAffineMatrix(Theta, [1;p;q]);
A.Sym = 1

A = A * [1 ; 1]
[ 2 1 ] * A

A.channels

end