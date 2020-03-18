function disp(N)
%%
%  File: disp.m
%  Directory: 1_PhD_projects/00_my_toolboxes/algo_P/@PGenAffineMatrix
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2018. October 21.
%
%%

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

    % fprintf('Variables of Pi: %s\n', stringify(N.pivars))

    fprintf('Size of the matrix:  %dx%d\n', N.q, N.m)
    
    if N.F_EVALUATED && isnumeric(N.channels_value)
        fprintf('Channels [value]:    %s \n', pcz_num2str(N.channels_value))
    end
    
    nr_channels = sprintf('[%d]',N.s);
    
    fprintf('Channels %6s:     %s (%s)\n', nr_channels, N.stringify(sym(N.channels)), class(N.channels))
    fprintf('Substitute into:     %s (%s)\n', N.stringify(N.vars), class(N.vars))
    fprintf('\n')
end