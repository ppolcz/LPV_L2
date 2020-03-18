function [Ni] = get_matrices(N,varargin)
%% get_matrices
%  
%  File: get_matrices.m
%  Directory: 1_PhD_projects/00_my_toolboxes/algo_P/@PGenAffineMatrix
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2018. October 22.
%

%%

if strcmp(N.type,'right')

    if isempty(varargin)
        new_channels = 1:N.s;
    else
        varargin = cellfun(@(v) { v(:).' }, varargin);
        new_channels = horzcat(varargin{:});
    end

    Ni = cell(1,numel(new_channels));
    
    for k = 1:numel(new_channels)
        Ni{k} = N.get_matrixk(k);
    end
    
else

    error('not implemented only for `right` typed matrices');

end
    
end