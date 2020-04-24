function varargout = subsref(N, S)
%% subsref
%  
%  File: subsref.m
%  Directory: 1_PhD_projects/00_my_toolboxes/algo_P/V11/@PAnnihilator
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2019.09.08. (szeptember  8, vasárnap), 13:07 [ATIRVA - OK]
%

%%
switch S(1).type
    case '()'

        % ([1;2],[3 4],5) --> {1 ; 2 ; 3 ; 4 ; 5}
        subs = cellfun(@(c) {c(:)},S.subs);
        x = vertcat(subs{:});
        
        if size(x,1) ~= numel(N.vars)
            display(x,pcz_OK_FAILED(0,'subsref input arguments'))
            error('You need to substitute into %s, but %d values provided (see above).', ...
            N.stringify(N.vars), size(x,1))
        end
        
        if strcmp('right',N.type)
            varargout{1} = N.Theta * kron(N.Im, N.xi(x));
        elseif strcmp('left',N.type)
            varargout{1} = N.Theta * kron(N.xi(x),N.Im);
        end
        
    case '{}'
        error('{} indexing not supported');
    case '.'
        % isprop(N,S(1).subs)
        if any(strcmp(fieldnames(N),S(1).subs))
            varargout{1} = builtin('subsref', N, S);
        elseif nargout > 0
            varargout = cell(1,nargout);
            [varargout{:}] = builtin('subsref', N, S);
        else
            builtin('subsref', N, S);
        end
end

end