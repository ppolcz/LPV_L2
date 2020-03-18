function [Nt] = transpose(N)
%% transpose
%  
%  File: transpose.m
%  Directory: 1_PhD_projects/00_my_toolboxes/algo_P/@PGenAffineMatrix
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  TODO
%

%%

    [~,Theta_cell] = N.get_matrices(N.channels);
    
    Theta_cell = cellfun(@transpose, Theta_cell, 'UniformOutput', 0);
    
    Nt = PGenAffineMatrix(Theta_cell, N.channels);

    Nt = Nt.copy(N);
    
end