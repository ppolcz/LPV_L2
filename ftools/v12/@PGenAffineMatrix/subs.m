function [M] = subs(N, vars, values, rvars)
%% subs
%  
%  File: subs.m
%  Directory: 7_ftools/ftools/v11/@PGenAffineMatrix
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2019. September 08. (2019a)
%

%%

if iscell(values)
    M = cellfun( ...
        @(v) { PGenAffineMatrix(N.Theta,subs(N.channels,vars,v),rvars) }, ...
        values );
else

    M = PGenAffineMatrix(N.Theta,subs(N.channels,vars,values),rvars);

end


end