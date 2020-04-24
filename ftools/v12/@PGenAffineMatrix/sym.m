function [ret] = sym(N)
%% sym
%  
%  File: sym.m
%  Directory: 1_PhD_projects/00_my_toolboxes/algo_P/@PGenAffineMatrix
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  TODO
%

%%

if ~N.issym
    N = N.generate_symbolic;
end

ret = N.Sym;

end