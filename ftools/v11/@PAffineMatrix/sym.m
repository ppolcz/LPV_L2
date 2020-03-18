function [ret] = sym(N)
%% sym
%  
%  File: sym.m
%  Directory: 1_PhD_projects/00_my_toolboxes/algo_P/@PAffineMatrix
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2018. October 22.
%

%%

if ~N.issym
    N = N.generate_symbolic;
end

ret = N.Sym;

end