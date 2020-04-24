function N = set_subsvars(N, subsvars)
%% set_subsvars
%  
%  File: set_subsvars.m
%  Directory: 1_PhD_projects/00_my_toolboxes/algo_P/@PAffineMatrix
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2018. October 22.
%

%%

    N.subsvars = subsvars;

    % disp =============================================
    % channels = N.channels
    % subsvars = N.subsvars
    % disp =============================================

    % if N.channels is an sdpvar
    symchannels = sym(N.channels);
    
    % 2019.02.14. (február 14, csütörtök), 16:29
    % [NEW METHOD]
    s1 = numel(symchannels);
    s2 = numel(subsvars);
    Act_per_New = repmat(symchannels(:),[1 s2]) ./ repmat(transpose(subsvars(:)),[s1,1]);
    Act_per_New(Act_per_New ~= 1) = 0;
    N.subsT = double(Act_per_New);
    % [NEW METHOD - END]

    % OLD METHOD
    % N.subsT = double(P_kiemel(symchannels, subsvars));
    
    symvar_of_Channels = symvar(symchannels);
    N.subsx0 = double(subs(symchannels, symvar_of_Channels, symvar_of_Channels*0));
    
    % disp ---------------------------------------------
    % subsT = N.subsT
    % disp ---------------------------------------------
    % symvar_of_Channels
    % subsx0 = N.subsx0
    % disp ---------------------------------------------

    assert(numel(N.subsx0) == size(N.subsT,1), ...
        sprintf('Size of subsx0 = %dx%d, Size of subsT = %dx%d',...
        size(N.subsx0), size(N.subsT)));
    
    assert(size(N.Theta,2) == N.m*size(N.subsx0,1),...
        sprintf('Size of Theta = %dx%d, Size of kron(I,subsx0) = %dx%d',...
        size(N.Theta), size(kron(N.Im,N.subsx0))));

    assert(size(N.subsT,2) == size(N.subsvars,1) || isempty(N.subsvars), ...
        sprintf('Size(subsT,2) = %d, size(N.subsvars,1) = %d', ...
        size(N.subsT,2), size(N.subsvars,1)));
                
    % N.subsind = pcz_select_symvar(subsvars, N.vars);   

end