function [Thetai,new_channels] = get_channels__(X, varargin)
%% get_channels__
%  
%  File: get_channels__.m
%  Directory: 1_PhD_projects/00_my_toolboxes/FinslerTools/v11/@PAffineMatrix
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2019. February 14.
%

%%

    % Generate selection matrix [BEGIN]

    varargin = cellfun(@(v) { transpose(v(:)) }, varargin);
    new_channels = transpose([ varargin{:} ]);

    zero = sym('__ZERO__','real');

    new_channels = sym(new_channels);
    new_channels(new_channels == 0) = zero;

    s1 = numel(X.channels);
    s2 = numel(new_channels);
    Act_per_New = repmat(X.channels(:),[1 s2]) ./ repmat(transpose(new_channels(:)),[s1,1]);
    Act_per_New(Act_per_New ~= 1) = 0;

    Selection_Matrix = double(Act_per_New);

    % Generate selection matrix [END]

    % Updated on 2019.11.05. (november  5, kedd), 21:15
    %
    %  Θ_new = Θ ( Im ⊗ E )
    %  E: Selection matrix
    %  eg. [ 0, 0]
    %      [ 1, 0]
    %      [ 0, 1]
    if isright(X)
        Thetai = X.Theta * kron(X.Im, Selection_Matrix);
    else
        Thetai = X.Theta * kron(Selection_Matrix, X.Im);
    end
        
end