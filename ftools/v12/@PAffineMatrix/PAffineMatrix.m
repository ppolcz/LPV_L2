classdef PAffineMatrix
%% PAffineMatrix
%
%  File: PAffineMatrix.m
%  Directory: 7_ftools/ftools/v11/@PAffineMatrix
%  Author: Peter Polcz (ppolcz@gmail.com) 
% 
%  Created on 2018. October 21.
%  Modified on 2019. November 05. (2019a): left, right, plfr
% 
% Affine matrix representation: 
% 
%  N(p) = [ N0 N1 ... Nnp ] * ([1;p] o Im) = ThetaR * (Im o [1;p]),
%         (right-typed representation)      (left-typed representation)
% 
% where `o' denotes the Kronecker tensor product and Im is an m x m
% identity matrix.
% 
% N0 and Ni corresponds to channels `1' and `pi', respectively, where p1,
% ..., pnp are free indefinit variables.
% 
% Dimensions of N(p): ny x nu (or q x m).
% 
% Nr. of pi parameters: np.
% 
% Nr. of channels: s = np + 1.


%%

properties (Constant = true)
    TYPE_RIGHT = 1;
    TYPE_LEFT = -1;
end

properties (GetAccess = public, SetAccess = private)

    % Variables appearing in the annihilator, eg. [x1,x2,p1,p2,dP1,dP2]
    vars

    % channels of the coefficient matrix, eg. [1,x1,x2,p1,p2,dP1,dP2]
    channels

    % Substitute into, eg. [x1,p1,p2]
    subsvars

    % By default each variable is set to 0 during the evaluation
    % (substitution). `subsx0' gives a zero vector, which corresponds to
    % vector `subsvars'.
    subsx0

    % It can happen, that the evaluation (substitution) is done in a
    % subspace of `vars' e.g. [x1,p1,p2] of vars [x1,x2,p1,p2,dP1,dP2].
    % In this case consider transformation matrix T:
    % 
    %  N.subsT = 
    %      0     0     0
    %      1     0     0
    %      0     0     0
    %      0     1     0
    %      0     0     1
    %      0     0     0
    %      0     0     0
    %
    % vars_VALUE = T * subsvars_VALUE.
    subsT

    symbolic = []

    Vertical_Sep = []
    Horizontal_Sep = []

    F_EVALUATED = 0;

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

    %{
    Coefficient matrix.
    Right-typed: N(p) = Theta * (I o [1;p])
    Left-typed: N(p) = Theta * ([1;p] o I), where Theta = [N0 N1 ... Nnp]
    %}
    Theta
end

properties (GetAccess = public, SetAccess = public)
    type = PAffineMatrix.TYPE_RIGHT;
    
    name = 'X';
    
    caption = '';

    Theta_right = [];
    Theta_left = [];
end

methods
    % function set.subsvars(N, subsvars_)
    
    % LEHET HOGY VEGTELEN CIKLUSBA KEVEREDNE :)
    % function N = set.channels(N,channels_)
    %     assert(numel(channels_) == numel(N.channels), ...
    %         'Nr. of new channels (%d), must coincide the nr. of old channels (%d)',...
    %         numel(N.channels), numel(channels_));
    %     N.channels = channels_;
    % end
    
    function Theta = get.Theta(N)
        if isright(N)
            Theta = N.Theta_right;
        else
            Theta = N.Theta_left;
        end
    end
        
    
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
        if ~issym(N)
            N = N.generate_symbolic;
        end
        Sym = N.symbolic;
    end
    
    function N = set.type(N,value)
        if ischar(value)
            if strcmp(value, 'right')
                N.type = PAffineMatrix.TYPE_RIGHT;
            elseif strcmp(value, 'left')
                N.type = PAffineMatrix.TYPE_LEFT;
            end
        elseif isscalar(value) && isnumeric(value)
            N.type = sign(value);
        end
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

        if isright(N)
            N.Theta_right = Theta_;
            N.Theta_left = [];
        else
            N.Theta_left = Theta_;
            N.Theta_right = [];
        end
        
        N = generate_other_Theta(N);
        
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
    
    function N = PAffineMatrix(varargin)
        % disp CONSTRUCTOR
        % display(nargin, 'Nr. of input arguments')
        
        % for i=1:nargin
        %     display(inputname(i), sprintf('inputname(%d)',i))
        % end
        
        N = N.init(varargin{:});
    end
    
    
    varargout = subsref(N, S)    
    [q,m] = size(N,dim);
    disp(N);
    disp_channels(N);
    
    N = double(N);
    N = value(N);
    
    % Mathematical operators
    ret = plus(X,Y);
    ret = mtimes(X,Y);
    ret = transpose(N);
    ret = ctranspose(N);
    ret = He(N);

    % Convert to a symbolic matrix
    ret = sym(N);
    sysplfr = plfr(N,pl_cell);
    syslfr = lfr(N,pl_cell);
    
    % Boolean operators:
    ret = isnan(N);

    [N_new] = to_symTheta(N,name,assumptions);
    
	N = set_subsvars(N, subsvars); % old
    N = set_vars(N, vars); % new (as in PGenAffineMatrix)
    N = set_channels(N, channels, vars);
    N = generate_symbolic(N);
    Ni = get_channels(N,varargin);
    [Ni,Ni_cell] = get_matrices(N,varargin);
    X = copy(X,Y)
            
    [X,Y] = adapt_channels(X,Y)

    function ret = issym(N)
        ret = ~isempty(N.symbolic);
    end
    
    function ret = stringify(N,s)
        ret = cell2mat(join(cellfun(@char, num2cell(s), 'UniformOutput', 0),', '));
    end
    
    function ret = isright(N)
        ret = N.type == PAffineMatrix.TYPE_RIGHT;
    end
        
    N = set_channel_by_index(N,Ni,i);

end

methods (Access = private)
    
    N = init(N, Theta_, channels_, subsvars_, varargin);
    
    %{
    If Theta is given as: { N0 N1 ... Nnp }.
    In this case Theta_right can be easily assigned as [ N0 N1 ... Nnp ]
    %}
	N = set_channels_from_cell(N,Theta_cell);
    [Thetai,new_channels] = get_channels__(N,varargin);
    
    %{
    Convert N(p) = Theta_right * (I o xi) = Theta_left * (xi o I),
    where Theta_left = [N0 N1 ... Nnp]
    %}
    N = generate_other_Theta(N);

end

end
