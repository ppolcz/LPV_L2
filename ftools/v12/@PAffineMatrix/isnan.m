function [ret] = isnan(N)
%% isnan
%  
%  File: isnan.m
%  Directory: 1_PhD_projects/00_my_toolboxes/FinslerTools/v11/@PAffineMatrix
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2019. February 14.
%

%%

ret = isdouble(N.Theta) && any(any(isnan(N.Theta)));

end