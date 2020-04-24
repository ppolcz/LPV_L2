function N = set_channel_by_index(N, Thetai, i)
%% set_channel_by_index
%  
%  File: set_channel_by_index.m
%  Directory: 1_PhD_projects/00_my_toolboxes/algo_P/@PAffineMatrix
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2018. October 22.
%

%%

    old_type = N.type;
    
    N.type = PAffineMatrix.TYPE_RIGHT;
    N.Theta_right(:,i:N.s:end) = Thetai;        
    N = N.generate_other_Theta;
    
    N.type = old_type;

end
