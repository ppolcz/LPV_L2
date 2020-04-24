function [Nt] = ctranspose(N)
%% transpose
%  
%  File: transpose.m
%  Directory: 1_PhD_projects/00_my_toolboxes/algo_P/@PAffineMatrix
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2018. October 22.
%

%%

    [~,Theta_cell] = N.get_matrices(N.channels);
    
    Theta_cell = cellfun(@ctranspose, Theta_cell, 'UniformOutput', 0);
    
    Nt = PAffineMatrix(Theta_cell, N.channels);
    
    Nt = Nt.copy(N);
    
end