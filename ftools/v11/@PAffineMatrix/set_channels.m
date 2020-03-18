function N = set_channels(N, channels)
%% set_channels
%  
%  File: set_channels.m
%  Directory: 1_PhD_projects/00_my_toolboxes/FinslerTools/v11/@PAffineMatrix
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2019. February 14.
%

%%

N = PAffineMatrix(N.get_channels__(channels),channels,N.subsvars);

end