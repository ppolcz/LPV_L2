function [ret] = mtimes(X,Y)
%% mtimes
%  
%  File: mtimes.m
%  Directory: 1_PhD_projects/00_my_toolboxes/algo_P/@PAffineMatrix
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2018. October 21.
%

%%

XisAffinMat = isa(X,'PAffineMatrix');
YisAffinMat = isa(Y,'PAffineMatrix');


if XisAffinMat && ~YisAffinMat
    
    ret_Theta = X.Theta * kron(Y,X.Is);
    
    ret = PAffineMatrix(ret_Theta,X.channels,X.subsvars);
    
    if X.issym, ret.Sym = 1; end
    
    return

elseif ~XisAffinMat && YisAffinMat
    
    ret_Theta = X * Y.Theta;
    
    ret = PAffineMatrix(ret_Theta,Y.channels,Y.subsvars);
    
    if Y.issym, ret.Sym = 1; end
    
    return
end


ret = [];

end
