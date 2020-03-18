function [r_indices,r_maxdiff,r_perc] = pcz_fhzero_report(z_fh, w, prec, N, varargin)
%% pcz_fhzero_report
%  
%  File: pcz_fhzero_report.m
%  Directory: 2_demonstrations/lib/matlab
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2018. June 09.

%%

if nargin < 4 || isempty(N) || ischar(N)
    if nargin >= 4 && ischar(N)
        varargin = [N varargin];
    end
    N = 10;
end

if nargin < 3 || isempty(prec) || ischar(prec)
    if nargin >= 3 && ischar(prec)
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

TMP_ZNEWEagSzRkCGbFsczkg =  pcz_dispFunctionName(title,'',struct('parent',1));

pcz_info('Tolerance: %g.', tol)

s = numel(w);

w_dummy = zeros(s,1);
z_dummy = z_fh(w_dummy);

p = numel(z_dummy);

ZERO = zeros(p,N);    
for i=1:N
    ZERO(:,i) = reshape(z_fh(rand(numel(w),1)), [p,1]);
end

greater = sum(abs(ZERO) > tol,2);
indices = find(greater);
perc = numel(indices) / p;
maxdiff = max(abs(ZERO(:)));

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

if nargout > 0
    r_indices = indices;
    r_maxdiff = maxdiff;
    r_perc = perc;
end

if p == 0
    perc = 0;
    maxdiff = 0;
end

if nargout == 0
    bool = perc == 0 && maxdiff < tol;
    
    pcz_dispFunction('Maximal difference: %g', maxdiff)
    pcz_dispFunctionSeparator
    
    pcz_info(bool, varargin{:}, {'first', first+1})
    
    if p == 0
        pcz_dispFunction('Variable is empty');
    end
    
    if ~bool
        pcz_dispFunction('Equality percentage: %g%%', (1-perc)*100)

        if ~isempty(indices)
            pcz_dispFunction('Indices, where not equal: %s / %d', pcz_num2str(indices(:)', 'format', '%d'), p);
            pcz_dispFunction('nr of elements where not zero/all: %d/%d', numel(indices), p)
        end

        % pcz_dispFunctionStackTrace
        % pcz_dispFunctionNeedNewLine
    end
end

pcz_dispFunctionEnd(TMP_ZNEWEagSzRkCGbFsczkg);

end