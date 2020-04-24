function [ret] = uminus(X)
%% uminus
%  
%  File: uminus.m
%  Directory: 1_PhD_projects/00_my_toolboxes/FinslerTools/v11/@PGenAffineMatrix
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  TODO
%

%%

ret = PGenAffineMatrix(-X.Theta,X.channels,X.subsvars);

if X.issym, ret.Sym = 1; end


end