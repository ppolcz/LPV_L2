function [r_indices,r_maxdiff,r_perc] = pcz_lfrzero_report(M_plfr, prec, N, varargin)
%% pcz_lfrzero_report
%  
%  File: pcz_lfrzero_report.m
%  Directory: utilities/symbolical
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2020. March 17. (2019b)
%

%%

if isa(M_plfr, 'lfr')
    M_plfr = plfr(M_plfr);
end
M_plfr = minlfr(M_plfr);

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

w_dummy = zeros(M_plfr.np,1);
z_dummy = M_plfr(w_dummy);

p = M_plfr.ny * M_plfr.nu;

if M_plfr.m1 > 0
    samples = generate(N, M_plfr.bounds, @(x) x);
    ZERO = reshape( M_plfr.val(samples) , [M_plfr.nu M_plfr.ny N] );
    ZERO_max = max(abs(ZERO),[],3);
else
    % 2020.04.10. (április 10, péntek), 17:46
    ZERO_max = abs(M_plfr.A);
end

% lfr2gss(N.lfr

% ZERO_mean = repmat(sum(ZERO,3) / N, [1 1 N]);
% ZERO_variance = sqrt(sum((ZERO - ZERO_mean).^2,3) / (N-1));

greater = ZERO_max(:) > tol;
indices = find(greater);
perc = numel(indices) / p;
maxdiff = max(ZERO_max(:));

S = dbstack;
% S.name

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
    
    pcz_dispFunction('Maximal error: %g', maxdiff)
    pcz_dispFunctionSeparator
    
    pcz_info(bool, varargin{:}, {'first', 2})
    
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


function samples = generate(N,lims,proj)

    if isempty(lims)
        samples = proj(randn(size(lims,1),N));
    else    
        samples_cell = cellfun(@(a,b) {rand(1,N)*(b-a)+a}, num2cell(lims(:,1)), num2cell(lims(:,2)));
        samples = proj(vertcat(samples_cell{:}));
        
        % After the projection some points may lie outside `lims'.
        samples = samples(:,prod(lims(:,1)*ones(1,N) <= samples & samples <= lims(:,2)*ones(1,N),1) == 1);
        
        %{
        plot3(samples(1,:),samples(2,:),samples(3,:),'.','MarkerSize',10)
        %}
        
    end

end
