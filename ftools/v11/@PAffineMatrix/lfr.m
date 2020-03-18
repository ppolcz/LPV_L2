function syslfr = lfr(N,pl_cell)
%% lfr
%  
%  File: lfr.m
%  Directory: 7_ftools/ftools/v11/@PAffineMatrix
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2019. November 05. (2019a)
%

%%

if nargin > 1 && iscell(pl_cell) && isa(pl_cell{1},'lfr')
    p = cellfun(@(p) { sym(p.blk.names{1}) }, pl_cell);
    p = vertcat(p{:});
    
    Xi_cell = cellfun(@(c) { N.Im * channel2lfr(c,p,pl_cell) }, num2cell(N.channels));    
else
    Xi_cell = cellfun(@(c) { N.Im * channel2lfr_simple(c) }, num2cell(N.channels));
end
    
Xi = vertcat(Xi_cell{:});

syslfr = N.Theta_left * Xi;

end

function c = channel2lfr(c,p,pl_cell)
    
    if isempty(symvar(c))
        c = lfr(double(c));
    elseif isscalar(symvar(c))
        c = pl_cell{p-c == 0};
    end
end

function c = channel2lfr_simple(c)
    if isempty(symvar(c))
        c = lfr(double(c));
    elseif isscalar(symvar(c))
        c = lfr(char(c));
    end
end