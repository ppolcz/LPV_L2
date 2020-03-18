function N = init(N, Theta_, channels_, vars_)
%% init
%  
%  File: init.m
%  Directory: 1_PhD_projects/00_my_toolboxes/algo_P/@PGenAffineMatrix
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2018. October 22.
%

%%

    % channels of the coefficient matrix, eg. [1,x1,x2,d1,x1*x2,x1*x3/(1+x1)]
    N.channels = channels_(:);

    % variables: [x1,x2,x3,d1]
    N.vars = vars_;

    if isa(Theta_, 'PGenAffineMatrix')
        
        N.Theta = Theta_.Theta;
        
    elseif iscell(Theta_)
    
        assert(numel(Theta_) == numel(channels_), ...
            sprintf(['If Theta is given as a cell, the number of matrices ' ...
            'in the cell must coincide the number of channels!\n' ...
            'Size of Theta:  %d x %d\n'...
            'Size of channels: %d x %d\n'], size(Theta_), size(channels_)));

        
        N.Theta = [ Theta_{:} ];

        N = N.set_matrices(Theta_);
    
    else

        N.Theta = Theta_;

    end
    
    % Substitute as..., eg. [x1,x1,d1,d2,Dd1,Dd2]
    N = N.set_vars(vars_);
    
end