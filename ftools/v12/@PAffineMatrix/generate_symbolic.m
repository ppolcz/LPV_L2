function N = generate_symbolic(N)
%% generate_symbolic
%  
%  File: generate_symbolic.m
%  Directory: 1_PhD_projects/00_my_toolboxes/algo_P/@PAffineMatrix
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2018. October 22.
%

%%

if isa(N.Theta,'double') || isa(N.Theta,'sym')
    
    if isright(N)
        N.symbolic = N.Theta * kron(N.Im,sym(N.channels));
    else
        N.symbolic = N.Theta * kron(sym(N.channels),N.Im);
    end
end

end