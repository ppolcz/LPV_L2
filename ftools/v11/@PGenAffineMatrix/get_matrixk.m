function [Ni] = get_matrixk(N,k)
%% get_matrices
%  
%  File: get_matrices.m
%  Directory: 1_PhD_projects/00_my_toolboxes/algo_P/@PGenAffineMatrix
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2018. October 22.
%

%%

if 0 < k && k <= numel(N.channels)

    if strcmp(N.type,'right')
    
        Ni = N.Theta(:,k:numel(N.channels):end);
        
    else
        
        Ni = 'Not impemented only for `right` typed matrices';
    
    end
    
else
    
    error('No channel nr. %d. Channels go from 1 to %d.', k, numel(N.channels))

end

end