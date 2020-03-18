function N = set_channels_from_cell(N, Theta_cell)
%% set_channels_from_cell
%  
%  File: set_channels_from_cell.m
%  Directory: 1_PhD_projects/00_my_toolboxes/algo_P/@PAffineMatrix
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2018. October 22.
%

%%

    assert(numel(Theta_cell) == numel(N.channels), ...
        'Set channels from cell of matrices: Number of matrices does not coincide to the number of channels');

    for i=1:numel(N.channels)
        N = N.set_channel_by_index(Theta_cell{i},i);
    end

end