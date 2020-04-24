function [getvalue] = G_VERBOSE(setvalue)
%% G_VERBOSE
%  
%  File: G_VERBOSE.m
%  Directory: 7_ftools/utilities/global
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2020. April 03. (2019b)
%

global VERBOSE

if nargin == 1
    VERBOSE = logical(setvalue);
elseif ~isscalar(VERBOSE)
    VERBOSE = false;
end

getvalue = VERBOSE;

end