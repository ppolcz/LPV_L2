function [r_indices,r_maxdiff,r_perc] = pcz_symzero_report(z, prec, N, varargin)
%% Script pcz_symzero_perc
%  
%  file:   pcz_symzero_perc.m
%  author: Peter Polcz <ppolcz@gmail.com> 
%  
%  Created on 2017.08.01. Tuesday, 01:21:16
%
%%

if nargin < 3 || isempty(N) || ischar(N)
    if nargin >= 3 && ischar(N)
        varargin = [N varargin];
    end
    N = 10;
end

if nargin < 2 || isempty(prec) || ischar(prec)
    if nargin >= 2 && ischar(prec)
        varargin = [prec varargin];
    end
    prec = 10;
end

if prec < 1
    tol = prec;
else
    tol = 10^(-prec);
end

title = '';
if ~isempty(varargin)
    title = sprintf(varargin{:});
end

s = symvar(sym(z));

if isempty(s)
    ZERO = double(z);
    indices = find(abs(ZERO(:)) > tol);
    perc = numel(indices) / numel(ZERO);
    maxdiff = max(abs(ZERO(:)));
    alldiff = abs(ZERO(:));
else
    z_fh = matlabFunction(z(:), 'vars', {s(:)});
    ZERO = zeros(numel(z),N);    
    for i=1:N
        ZERO(:,i) = z_fh(rand(numel(s),1));
    end
    
    greater = sum(abs(ZERO) > tol,2);
    indices = find(greater);
    perc = numel(indices) / numel(z);
    maxdiff = max(abs(ZERO(:)));
    alldiff = abs(ZERO(:));
end

S = dbstack;
% S.name

first = 1;
for i = 1:numel(S)-1
    if strcmp(S(2).name, 'pcz_symzero_report')
        first = i;
        break;
    end
end

if numel(S) > first && strcmp(S(first+1).name, 'pcz_symeq_report')
    first = first+1;
end

if nargout > 2
    r_perc = perc;
end

if nargout > 1
    r_maxdiff = maxdiff;
end

if nargout > 0
    r_indices = indices;
end

if nargout == 0
    TMP_ZNEWEagSzRkCGbFsczkg =  pcz_dispFunctionName(title,'',struct('parent',1));
    pcz_info('Tolerance: %g.', tol)

    bool = perc == 0 && maxdiff < tol;
    pcz_info('Maximal difference: %g', maxdiff);
    
    pcz_dispFunctionSeparator
    pcz_info(bool, varargin{:}, {'first', first+1})
        
    if ~bool
        pcz_dispFunction('Equality percentage: %g%%', (1-perc)*100)

        if ~isempty(indices)
            pcz_dispFunction('Indices, where not equal: %s / %d', pcz_num2str(indices(:)','format', '%d'), numel(z));
            pcz_dispFunction('nr of elements where not zero/all: %d/%d', numel(indices), numel(z));
            pcz_dispFunction('differences: %s', pcz_num2str(alldiff(indices), 'format', '%g'))
        end

        % pcz_dispFunctionStackTrace
        % pcz_dispFunctionNeedNewLine
    end
    
    pcz_dispFunctionEnd(TMP_ZNEWEagSzRkCGbFsczkg);
end
%%

