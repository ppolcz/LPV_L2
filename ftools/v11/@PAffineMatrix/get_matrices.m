function [Ni,Ni_cell] = get_matrices(N,varargin)
%% get_matrices
%  
%  File: get_matrices.m
%  Directory: 1_PhD_projects/00_my_toolboxes/algo_P/@PAffineMatrix
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2018. October 22.
%  Modified on 2019. November 05. (2019a)
%

%%

    if nargin == 1
        % Simpler: 2019.11.05. (november  5, kedd), 21:24

        %{
        N = PAffineMatrix([1 2 3 4 5 6],[1;sym('x')])

        Assume that m = 3 (cols of N(p)), s = 2 (nr. of channels)

        Theta_right = 
            p1 p2 | p1 p2 | p1 p2
            1  2  | 3  4  | 5  6 

        Theta_left = 
            p1    | p2
            1 3 5 | 2 4 6

        Then, I =
             1     2
             3     4
             5     6
        %}
        I = (1:N.s:N.s*N.m)'+(0:N.s-1);

        I_cell = num2cell(I,1);
        Ni_cell = cellfun(@(I) {N.Theta(:,I)}, I_cell);
    
    else
        % Older, but more general
    
        varargin = cellfun(@(v) { transpose(v(:)) }, varargin);
        new_vars = [ varargin{:} ];
    
        % Number of channels in `Thetai` coefficient matrix with respected
        % to vars
        sp = numel(new_vars);

        N_channels = N.get_channels(new_vars);

        Ni_cell = cell(1,sp);

        for i = 1:sp
            Ni_cell{i} = N_channels.Theta(:,i:sp:end);
        end
        
    end

    Ni = [Ni_cell{:}];

end