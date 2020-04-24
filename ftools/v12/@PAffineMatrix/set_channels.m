function N = set_channels(N, channels_, vars_)
%% set_channels
%  
%  File: set_channels.m
%  Directory: 1_PhD_projects/00_my_toolboxes/FinslerTools/v11/@PAffineMatrix
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2019. February 14.
%

%%

if nargin == 2
    
    % Old operation mode (for compatibility reasonse)
    N = PAffineMatrix(N.get_channels__(channels_),channels_,N.subsvars);

elseif nargin == 3
    
    % New operation mode (same as in PGenAffineMatrix)

    if N.s ~= numel(channels_)
        error('Nr. of channels (%d) must conincide to the nr. of requested channels (%d)',...
            N.s, numel(channels_));
    end

    N.channels = channels_;
    N = N.set_vars(vars_);
end
    
end