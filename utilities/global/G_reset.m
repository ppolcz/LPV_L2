function G_reset(verbose)
%% pcz_reset
%  
%  File: pcz_reset.m
%  Directory: 7_ftools/utilities/global
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2020. April 03. (2019b)
%

global SCOPE_DEPTH LATEX_EQNR 

if nargin < 1
    verbose = true;
end

G_VERBOSE(verbose);

SCOPE_DEPTH = 0;
LATEX_EQNR = 0;

end