function [ret] = uminus(X)
%% uminus
%  
%  File: uminus.m
%  Directory: 1_PhD_projects/00_my_toolboxes/FinslerTools/v11/@PAffineMatrix
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2019. February 14.
%

%%

ret = PAffineMatrix(-X.Theta,X.channels,X.subsvars);

if X.issym, ret.Sym = 1; end


end