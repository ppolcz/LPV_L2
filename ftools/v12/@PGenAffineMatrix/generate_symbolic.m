function N = generate_symbolic(N)
%% generate_symbolic
%  
%  File: generate_symbolic.m
%  Directory: 1_PhD_projects/00_my_toolboxes/algo_P/@PGenAffineMatrix
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2018. October 22.
%

%%

if isa(N.Theta,'double') || isa(N.Theta,'sym')
    N.symbolic = N.Theta * kron(N.Im,sym(N.channels(:)));
end

end