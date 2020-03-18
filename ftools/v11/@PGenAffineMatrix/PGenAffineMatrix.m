classdef PGenAffineMatrix
%% PGenAffineMatrix
%
%  File: PGenAffineMatrix
%  Directory: 1_PhD_projects/00_my_toolboxes/algo_P/@PGenAffineMatrix
%  Author: Peter Polcz (ppolcz@gmail.com) 
% 
%  Created on 2018. October 21.


%%

properties (Constant = true)
end

properties (GetAccess = public, SetAccess = private)

    % Variables [x1,x2,p1,p2,dp1,dp2,...]
    vars
    
    % Channels of the coefficient matrix, eg. [1,x1,x1*x2,p1*p2*x1,sin(x1),...]
    channels

    symbolic = []

    Vertical_Sep = []
    Horizontal_Sep = []
    
    F_EVALUATED = 0;
    
    % Coefficient matrix
    Theta = []
    xi
end

properties (Dependent)
    Header
    Table
    Matrix
    Sym
    
    s % number of channels
    m,nu % number of columns of the matrix
    q,ny % number of rows of the matrix

    Is,Im
    
    channels_value
end

properties (GetAccess = public, SetAccess = public)
    type = 'right';
    name = 'X'
    caption = '';
end

methods
    
    function Header = get.Header(N)
        c = sym('I','real');
        Header = repmat(transpose([ c ; N.vars ]),[1 N.m]);
    end
    
    function Table = get.Table(N)
        Table = [ N.Header ; sym(N.Theta) ];
        if ~N.issym
            Table = vpa(Table,3);
        end
    end
    
    function Matrix = get.Matrix(N)
        Matrix = N.Theta * kron(N.Im, N.channels);
    end
    
    function Sym = get.Sym(N)
        Sym = N.symbolic;
    end
    
    function N = set.Sym(N,value)
        if value == 1
            N = N.generate_symbolic;
        elseif isempty(value) || value == 0
            N.symbolic = [];
        end
    end
    
    function N = set.Theta(N,Theta_)
        assert(isempty(N.Theta) || all(size(Theta_) == size(N.Theta)), ...
            sprintf('Size(old_Theta) = %dx%d != %dx%d = Size(new_Theta)',...
            size(N.Theta), size(Theta_)));
        
        N.Theta = Theta_;
        
        if N.issym
            N = N.generate_symbolic;
        end
    end
    
    function s = get.s(N)
        s = numel(N.channels);
    end
    
    function m = get.m(N)
        m = size(N.Theta,2) / numel(N.channels);
    end
    
    function nu = get.nu(N)
        nu = N.m;
    end
    
    function q = get.q(N)
        q = size(N.Theta,1);
    end
    
    function ny = get.ny(N)
        ny = N.q;
    end
    
    function Is = get.Is(N)
        Is = eye(N.s);
    end
    
    function Im = get.Im(N)
        Im = eye(N.m);
    end
    
    function channels_value = get.channels_value(N)
        channels_value = value(N.channels);
    end
end

methods (Access = public)
    
    function N = PGenAffineMatrix(varargin)
        N = N.init(varargin{:});
    end
    
    
    varargout = subsref(N, S)    
    [q,m] = size(N,dim);
    disp(N);
    disp_channels(N);
    
    N = double(N);
    N = value(N);
    
    ret = plus(X,Y);
    ret = mtimes(X,Y);
    ret = sym(N);
    ret = transpose(N);
    ret = ctranspose(N);
    ret = isnan(N);
    ret = He(N);
    
    [N_new] = to_symTheta(N,name,assumptions);
    
    % Use it for partial substitution
    % vars: variables, which are ment to be substituted by their actual vaue
    % values: values
    % rvars: remaining values to substitute
    N = subs(N,vars,values,rvars);
    
	N = set_vars(N, vars);
    N = set_matrices(N, Theta_cell,channels);
    N = set_matrixk(N,Thetak,k);
    N = generate_symbolic(N);
    Theta = get_channels(N,varargin);
    Ni = get_matrices(N,varargin);
    Ni = get_matrixk(N,k);
    X = copy(X,Y)
            
    [X,Y] = adapt_channels(X,Y)

    function ret = issym(N)
        ret = ~isempty(N.symbolic);
    end
    
    function ret = stringify(N,s)
        ret = cell2mat(join(cellfun(@char, num2cell(s), 'UniformOutput', 0),', '));
    end
    
end

methods (Access = private)
    
    N = init(N, Theta_, channels_, vars_);
    
end

end
