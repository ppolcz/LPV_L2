function [r_indices,r_maxdiff,r_perc] = pcz_gsszero_report(M_gss, prec, N, varargin)
%% pcz_lfrzero_report
%  
%  File: pcz_lfrzero_report.m
%  Directory: utilities/symbolical
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2020. March 17. (2019b)
%

%%

if isa(M_gss,'plfr')
    M_gss = M_gss.lfrtbx_obj;
end

if isa(M_gss,'lfr')
    try
        M_gss = lfr2gss(M_gss);
    catch e
        if size(M_gss,'order') == 0
            M_gss = gss(M_gss.d);
        else
            error(e)
        end
    end          
end

if nargin < 3 || isempty(N) || ischar(N)
    if nargin >= 3 && ( ischar(N) || ischar(prec) )
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

TMP_ZNEWEagSzRkCGbFsczkg =  pcz_dispFunctionName(title,'',struct('parent',1));

pcz_info('Tolerance: %g.', tol)

if size(M_gss,'blk') == 0
    ZERO_max = M_gss.M.d;
else
    ZERO = dbsample(M_gss,N);
    ZERO = cat(3,ZERO{:});
    ZERO_max = max(abs(ZERO),[],3);
end

% lfr2gss(N.lfr

% ZERO_mean = repmat(sum(ZERO,3) / N, [1 1 N]);
% ZERO_variance = sqrt(sum((ZERO - ZERO_mean).^2,3) / (N-1));

[ny,nu] = size(M_gss);

greater = ZERO_max(:) > tol;
indices = find(greater);
perc = numel(indices) / (ny*nu);
maxdiff = max(ZERO_max(:));

S = dbstack;
% S.name

if nargout > 0
    r_indices = indices;
    r_maxdiff = maxdiff;
    r_perc = perc;
end

if ny*nu == 0
    perc = 0;
    maxdiff = 0;
end

if nargout == 0
    bool = perc == 0 && maxdiff < tol;
    
    pcz_dispFunction('Maximal difference: %g', maxdiff)
    pcz_dispFunctionSeparator
    
    pcz_info(bool, varargin{:}, {'first', 2})
    
    if ny*nu == 0
        pcz_dispFunction('Variable is empty');
    end
    
    if ~bool
        pcz_dispFunction('Equality percentage: %g%%', (1-perc)*100)

        if ~isempty(indices)
            pcz_dispFunction('Indices, where not equal: %s / %d', pcz_num2str(indices(:)', 'format', '%d'), ny*nu);
            pcz_dispFunction('nr of elements where not zero/all: %d/%d', numel(indices), ny*nu)
        end

        % pcz_dispFunctionStackTrace
        % pcz_dispFunctionNeedNewLine
    end
end

pcz_dispFunctionEnd(TMP_ZNEWEagSzRkCGbFsczkg);

end