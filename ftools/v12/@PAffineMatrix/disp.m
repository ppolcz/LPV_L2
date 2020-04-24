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
        if numel(N.Theta) > 200
            fprintf('Theta: [%dx%d], values between (%g, %g)\n', size(N.Theta), min(N.Theta(:)), max(N.Theta(:)))
        else
            % 2020.03.31. (március 31, kedd), 00:24
            if N.isright
                fprintf('Theta: [%dx%d] = [(%dx%d)x%d]\n', N.ny, N.nu*N.s, N.ny, N.s, N.nu);
                indices = reshape(1:N.s*N.m,[N.s,N.m]);
            else
                fprintf('Theta: [%dx%d] = [(%dx%d)x%d]\n', N.ny, N.nu*N.s, N.ny, N.nu, N.s);
                indices = reshape(1:N.s*N.m,[N.m,N.s]);
            end
            indices_cell = num2cell(indices,1);
            blocks_cell = cellfun(@(ind) {N.Theta(:,ind)}, indices_cell);

            blocks_str_cell = cellfun(@(block) { split(evalc('disp(block)'),newline) }, blocks_cell);
            blocks_str_cell = num2cell(horzcat(blocks_str_cell{:}),2);

            separator = ' ,';
            disp_line_by_line = cellfun(@(c) { strjoin(c,separator) }, blocks_str_cell);

            % get rid of empty rows:
            row_length = cellfun(@numel, disp_line_by_line);

            disp_line_by_line = disp_line_by_line(row_length > numel(separator) * numel(blocks_cell));

            disp_all = strjoin(disp_line_by_line, newline);

            disp(disp_all);
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