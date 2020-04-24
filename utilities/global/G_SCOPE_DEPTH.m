function [ret] = G_SCOPE_DEPTH(incvalue)
%% G_SCOPE_DEPTH
%  
%  File: G_SCOPE_DEPTH.m
%  Directory: 7_ftools/utilities/global
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2020. April 03. (2019b)
%

%%

global SCOPE_DEPTH

if ~isscalar(SCOPE_DEPTH)
    SCOPE_DEPTH = 0;
end

if nargin < 1
    ret = SCOPE_DEPTH;
    return
end

if incvalue == 0
    SCOPE_DEPTH = 0;
elseif abs(incvalue) == 1
    SCOPE_DEPTH = SCOPE_DEPTH + incvalue;
end
    

if SCOPE_DEPTH < 0
    warning 'SCOPE_DEPTH is negative, setting to 0'
    SCOPE_DEPTH = 0;
end


[ST,I] = dbstack;

if SCOPE_DEPTH - numel(ST) + I - 1 > 10
    warning 'SCOPE_DEPTH decreased'
    SCOPE_DEPTH = numel(ST) - I + 1;
end

ret = SCOPE_DEPTH;

end