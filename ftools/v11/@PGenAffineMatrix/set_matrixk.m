function N = set_matrixk(N, Thetak, k)
%% set_channel_by_index
%  
%  File: set_channel_by_index.m
%  Directory: 1_PhD_projects/00_my_toolboxes/algo_P/@PGenAffineMatrix
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2019.09.08. (szeptember  8, vasárnap), 13:07 [ATIRVA - OK]
%

%%

    nr = numel(N.channels);

    N.Theta(:,k:nr:end) = Thetak;        

end