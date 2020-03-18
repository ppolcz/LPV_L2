function [ret] = mtimes(X, Y)
%% mtimes
%  
%  File: mtimes.m
%  Directory: 1_PhD_projects/00_my_toolboxes/FinslerTools/v11/@plfr
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2019. April 10.
%

%%

X_var = [];
if isa(X,'plfr')
    X_var = X.subsvars(:);
    X = X.lfrtbx_obj;
elseif isa(X,'lfr')
    X_plfr = plfr(X);
    X_var = X_plfr.subsvars(:);
end


Y_var = [];
if isa(Y,'plfr')
    Y_var = Y.subsvars(:);
    Y = Y.lfrtbx_obj;
elseif isa(Y,'lfr')
    Y_plfr = plfr(Y);
    Y_var = Y_plfr.subsvars(:);
end

ret = plfr(X*Y,plfr.var_helper__3([X_var ; Y_var]));

end
