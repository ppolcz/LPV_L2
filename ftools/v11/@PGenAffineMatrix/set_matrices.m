function N = set_matrices(N, Theta_cell, channels_)
%% set_channels_from_cell
%  
%  File: set_channels_from_cell.m
%  Directory: 1_PhD_projects/00_my_toolboxes/algo_P/@PGenAffineMatrix
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2019.09.08. (szeptember  8, vasárnap), 13:07 [ATIRVA - OK]
%

%%

    if nargin < 3
        channels_ = 1:numel(N.channels);
    end

    if numel(Theta_cell) ~= numel(N.channels)
        error('Set channels from cell of matrices: Number of matrices (%d) does not coincide to the number of channels (%d)',...
            numel(Theta_cell), numel(N.channels));
    end

    for k = channels_
        N = N.set_matrixk(Theta_cell{k},k);
    end

end