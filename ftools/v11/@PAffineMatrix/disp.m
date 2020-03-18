function disp(N)
%%
%  File: disp.m
%  Directory: 1_PhD_projects/00_my_toolboxes/algo_P/@PAffineMatrix
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2018. October 21.
%
%%

    % Hozzaadva: 2019.11.05. (november  5, kedd), 21:06
        types = {'LEFT','RIGHT'};
        defin = {
            'Theta*(xi o Im), where Theta = [N1,...,Ns], dim(xi) = '
            'Theta*(Im o xi), where dim(xi) = '
            };
        typeind = (N.type+1)/2+1;

        fprintf('%s object, %s-typed: %s = %s%d\n',...
            class(N),types{typeind},N.name,defin{typeind},N.s);
    % Hozzaadás vége
    

    if N.issym
        fprintf('Symbolic form of the matrix: [%dx%d] \n', size(N))
        disp(N.Sym);
    elseif isa(N.Theta, 'sdpvar')
        fprintf('Theta: ');
        display(N.Theta);
    else
        fprintf('Theta: [%dx%d]', size(N.Theta));
        if numel(N.Theta) > 200
            fprintf(', values between (%g, %g)\n', min(N.Theta(:)), max(N.Theta(:)))
        else
            disp ' '
            disp(N.Theta);
        end
    end

    fprintf('Size of the matrix:  %dx%d\n', N.q, N.m)
    
    % Kikommentezve: 2019.11.05. (november  5, kedd), 21:06
    % if N.F_EVALUATED && isnumeric(N.channels_value)
    %     fprintf('Channels [value]:    %s \n', pcz_num2str(N.channels_value))
    % end
    
    fprintf('Channels (xi):       %s (%s)\n', N.stringify(sym(N.channels)), class(N.channels))
    fprintf('Substitute into:     %s (%s)\n', N.stringify(N.subsvars), class(N.subsvars))
    fprintf('\n')
end