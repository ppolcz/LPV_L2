function N = init(N, Theta_, channels_, subsvars_, varargin)
%% init
%  
%  File: init.m
%  Directory: 1_PhD_projects/00_my_toolboxes/algo_P/@PAffineMatrix
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2018. October 22.
%
% Before the assignment of Theta it is very important to set `type',
% because the assignment in `set.Theta' depends on type = left/right.


%% Parse arguments passed to init

    % In the minimal case only Theta_ and channels_ are given
    narginchk(3,Inf);
        
    % If subsvars_ is not given, but optional arguments are given
    if nargin > 3 && ischar(subsvars_)
        varargin = [subsvars_ varargin];
        subsvars_ = [];
    end
    
    N = parsepropval(N,varargin{:});
    
    
%% Set type and Theta
    
    % channels of the coefficient matrix, eg. [1,x1,x2,d1]
    N.channels = channels_;

    if isa(Theta_, 'PAffineMatrix')
        
        % Inherit type: left/right
        N.type = Theta_.type;
        
        % This will call set.Theta, where the assignment depends on type.
        N.Theta = Theta_.Theta;
        
    elseif iscell(Theta_)
    
        % Theta is given as { N0 N1 ... Nnp }
        assert(numel(Theta_) == numel(channels_), ...
            sprintf(['If Theta is given as a cell, the number of matrices ' ...
            'in the cell must coincide the number of channels!\n' ...
            'Size of Theta:  %d x %d\n'...
            'Size of channels: %d x %d\n'], size(Theta_), size(channels_)));

        % Store type given as an argument.
        requested_type = N.type;

        % Set type to LEFT
        N.type = PAffineMatrix.TYPE_LEFT;
        
        % Call set.Theta with a LEFT-typed Theta value
        N.Theta = [ Theta_{:} ];

        % Reset type as originally requested
        N.type = requested_type;
        
    else

        % In this case type is set as it was given in the arguments or by
        % its default value.
        
        % Call set.Theta
        N.Theta = Theta_;

    end        

%% Set channels and other symbolic variables
    
    % Variables appearing in the matrix, eg. [x1,x2,d1]
    channels_ = sym(channels_);
    channels_(channels_ == 1 | channels_ == 0) = [];
    N.vars = channels_;

    if nargin <= 3 || isempty(subsvars_)
        subsvars_ = N.vars;
    end

    % Substitute as..., eg. [x1,x1,d1,d2,Dd1,Dd2]
    N = N.set_subsvars(sym(subsvars_));
    
end