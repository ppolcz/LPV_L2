function [Nt] = ctranspose(N)
%% transpose
%  
%  File: transpose.m
%  Directory: 1_PhD_projects/00_my_toolboxes/algo_P/@PGenAffineMatrix
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2018. October 22.
%

%%

    Theta_cell = N.get_matrices();
    
    Theta_cell = cellfun(@ctranspose, Theta_cell, 'UniformOutput', 0);
    
    Nt = PGenAffineMatrix(Theta_cell, N.channels, N.vars);
    
    Nt = Nt.copy(N);
    
end