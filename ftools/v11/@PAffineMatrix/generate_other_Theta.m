function N = generate_other_Theta(N)
%% generate_other_Theta
%  
%  File: generate_other_Theta.m
%  Directory: 7_ftools/ftools/v11/@PAffineMatrix
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2019. November 05. (2019a)
%

%%
        
    if isright(N)

        assert(~isempty(N.Theta_right), ...
            'PAffineAnnihilator (name: `%s'') is RIGHT-typed, but Theta_right is empty!', ...
            N.name)

        assert(size(N.Theta_right,2) == N.s*N.m, ...
            'PAffineAnnihilator (name: `%s'', dim: %dx%d, nr. channels: %d) is RIGHT-typed, but size(Theta_right) = %dx%d', ...
            N.name, N.ny, N.nu, N.s, size(N.Theta_right))

        %{
        Assume that m = 3 (cols of N(p)), s = 2 (nr. of channels)

        Theta_right = 
            p1 p2 | p1 p2 | p1 p2
            1  2  | 3  4  | 5  6 

        Theta_left = 
            p1    | p2
            1 3 5 | 2 4 6

        Then, I =
             1     2
             3     4
             5     6
        %}

        I = (1:N.s:N.s*N.m)'+(0:N.s-1);
        N.Theta_left = N.Theta_right(:,I(:));

    else

        assert(~isempty(N.Theta_left), ...
            'PAffineAnnihilator (name: `%s'') is LEFT-typed, but Theta_left is empty!', ...
            N.name)

        assert(size(N.Theta_left,2) == N.s*N.m, ...
            'PAffineAnnihilator (name: `%s'', dim: %dx%d, nr. channels: %d) is LEFT-typed, but size(Theta_left) = %dx%d', ...
            N.name, N.ny, N.nu, N.s, size(N.Theta_left))

        %{
        Assume that m = 3 (cols of N(p)), s = 2 (nr. of channels)

        Theta_left = 
            p1    | p2
            1 2 3 | 4 5 6

        Theta_right = 
            p1 p2 | p1 p2 | p1 p2
            1  4  | 2  5  | 3  6 

        Then, I =
             1     2     3
             4     5     6
        %}

        I = (1:N.m:N.s*N.m)'+(0:N.m-1);
        N.Theta_right = N.Theta_left(:,I(:));

    end
end