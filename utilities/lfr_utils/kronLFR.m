function [XI] = kronLFR(xi,I)
%% kron
%  
%  File: kron.m
%  Directory: utilities/lfr_utils
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2019. June 05. (2019a)
%

%%

if isa(xi,'cell')

    XI_cell = cellfun( @(var) {var * I}, xi );
    
    XI = vertcat(XI_cell{:});
    
elseif isa(xi,'lfr')
    
    error 'Not implemented yet';
        
end

error 'Not supported';
        

end