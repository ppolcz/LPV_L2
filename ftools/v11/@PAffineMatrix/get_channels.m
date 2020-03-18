function Xi = get_channels(X,varargin)
%%
%  File: get_channels
%  Directory: /1_PhD_projects/00_my_toolboxes/algo_P/@PAffineMatrix
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2018. October 22.
% 
%  First: 2018.10.03. (október  3, szerda), 00:10
%
%  Θ_new = Θ ( Im ⊗ E )
%  E: Selection matrix
%  eg. [ 0, 0]
%      [ 1, 0]
%      [ 0, 1]
% 
% Usage: 
% 
% Nb.get_channels(1,d,x)
% 
%%

    [Thetai,new_channels] = X.get_channels__(varargin{:});
    
    % Updated on 2019.11.05. (november  5, kedd), 21:15
    Xi = PAffineMatrix(Thetai,new_channels,X.subsvars,'type',X.type);

end
