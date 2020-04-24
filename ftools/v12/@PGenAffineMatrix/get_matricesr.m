function [Thetai,new_channels] = get_matricesr(N, varargin)
%% get_channels__
%  
%  File: get_channels__.m
%  Directory: 1_PhD_projects/00_my_toolboxes/FinslerTools/v11/@PGenAffineMatrix
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2019.09.08. (szeptember  8, vasárnap), 13:07 [ATIRVA - OK]
% 
%  Selection matrix = 
%  eg. [ 0, 0]
%      [ 1, 0]
%      [ 0, 1]
% 

%%


    if strcmp(N.type,'right')
    
        varargin = cellfun(@(v) { v(:).' }, varargin);
        new_channels = horzcat(varargin{:});

        Mch1 = repmat((1:N.s)',[1 numel(new_channels)]);
        Mch2 = repmat(new_channels,[N.s 1]);

        Selection_Matrix = double(Mch1 == Mch2);

        Thetai = N.Theta * kron(N.Im, Selection_Matrix);
    
    else
    
        Thetai = 'not implemented yet';
        
    end
    
end