function N = set_channels(N, channels_, vars_)
%% set_channels_from_cell
%  
%  File: set_channels_from_cell.m
%  Directory: 1_PhD_projects/00_my_toolboxes/algo_P/@PGenAffineMatrix
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2019.09.08. (szeptember  8, vasárnap), 13:07 [ATIRVA - OK]
%

%%

if N.s ~= numel(channels_)
    error('Nr. of channels (%d) must conincide to the nr. of requested channels (%d)',...
        N.s, numel(channels_));
end

N.channels = channels_;
N = N.set_vars(vars_);

end