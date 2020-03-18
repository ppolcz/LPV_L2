function [varargout] = pcz_split_matrix(M, sizes_i, sizes_j, varargin)
%% pcz_split_matrix
%  
%  File: pcz_split_matrix.m
%  Directory: 2_demonstrations/lib/matlab
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2018. November 09.
%

opts.RowWise = true;
opts.ReturnCell = false;
opts = parsepropval(opts,varargin{:});

%%

% warning 'pcz_split_matrix: FIGYELEM változás történt! Régi működés: ''RowWise'', false'
% pcz_dispFunctionStackTrace('first',0,'last',0)

if isempty(sizes_i)
    sizes_i = size(M,1);
end

if isempty(sizes_j)
    sizes_j = size(M,2);
end

n = numel(sizes_i);
m = numel(sizes_j);

if (nargout > n*m)
    error('To many output arguments. Must not exceed %d x %d = %d', n, m, n*m);
end

% [1] HANDLE Inf in sizes: [3 2 Inf 4]

    inf_i = isinf(sizes_i);
    inf_j = isinf(sizes_j);

    % [3 2 Inf 4] es size(M) = [20,30] --> [3 2 *11* 4]
    sizes_i(inf_i) = ( size(M,1) - sum(sizes_i(~inf_i)) ) / sum(inf_i);
    sizes_j(inf_j) = ( size(M,2) - sum(sizes_j(~inf_j)) ) / sum(inf_j);

% [1:END]


ind_i = num2cell([0 cumsum(sizes_i)]);
ind_j = num2cell([0 cumsum(sizes_j)]);

to_domains = @(sizes) cellfun( @(from,to) {from+1:to}, ...
    sizes(1:end-1), sizes(2:end));

range_i = to_domains(ind_i);
range_j = to_domains(ind_j);

varargout = cell(n,m);

if opts.RowWise

    % Uj mukodese
    % Output arguments row-vise:
    % [A,B,C,D] = ~([A B ; C D], [nx nu], [nx ny])
    for i = 1:n
        for j = 1:m
            varargout{j+m*(i-1)} = M(range_i{i}, range_j{j});
        end
    end
    
else
    
    % Regi mukodese:
    % Output arguments column-vise:
    % [A,C,B,D] = ~([A B ; C D], [nx nu], [nx ny])
    for i = 1:n
        for j = 1:m
            varargout{i+n*(j-1)} = M(range_i{i}, range_j{j});
        end
    end
    
end

end


function test1
%%
    global SCOPE_DEPTH VERBOSE LATEX_EQNR 
    SCOPE_DEPTH = 0;
    VERBOSE = 1;
    LATEX_EQNR = 0;

    eq = @(A,B) all(all(A == B));
    
    [A,B,C,D] = pcz_split_matrix([1 2 ; 3 4], [1 1], [1 1]);
    pcz_OK_FAILED(A == 1 && B == 2 && C == 3 && D == 4, 'simple test\n')
    
    [A,C,B,D] = pcz_split_matrix([1 2 ; 3 4], [1 1], [1 1], 'RowWise', false);
    pcz_OK_FAILED(A == 1 && B == 2 && C == 3 && D == 4, 'simple test RowWise\n')

    A = pcz_split_matrix([1 2 3 ; 4 5 6 ; 7 8 9], [2], [2]);
    pcz_OK_FAILED(eq(A,[1 2 ; 4 5]), 'simple test, single matrix\n')
    
    [A,B,C,D] = pcz_split_matrix([1 2 3 ; 4 5 6 ; 7 8 9], [2 Inf], [2 Inf]);
    pcz_OK_FAILED(eq(A,[1 2 ; 4 5]) && eq(B,[3 ; 6]) && eq(C,[7 8]) && D==9, ...
        'simple test, with Inf\n')
    
    [A,B,C,D] = pcz_split_matrix([1 2 ; 3 4], [Inf Inf], [1 Inf]);
    pcz_OK_FAILED(A == 1 && B == 2 && C == 3 && D == 4, 'simple test, multiple Inf\n')
    
end
