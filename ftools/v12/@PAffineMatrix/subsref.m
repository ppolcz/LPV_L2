function varargout = subsref(N, S)
%% subsref
%  
%  File: subsref.m
%  Directory: 1_PhD_projects/00_my_toolboxes/algo_P/V11/@PAnnihilator
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2018. June 10.
%

%%

switch S(1).type
    case '()'
        
        % [REGI] 2018. June 10.
        % x = S.subs{1};
        
        % [UJ] 2019.11.05. (november  5, kedd), 19:18
        % ([1;2],[3 4],5) --> {1 ; 2 ; 3 ; 4 ; 5}
        subs = cellfun(@(c) {c(:)},S.subs);
        x = vertcat(subs{:});

        assert(size(x,1) == numel(N.subsvars), sprintf('You need to substitute into %s', N.stringify(N.subsvars)))
        
        if isright(N)
            varargout{1} = N.Theta * kron(N.Im, N.subsx0 + N.subsT*x);
        else
            varargout{1} = N.Theta * kron(N.subsx0 + N.subsT*x,N.Im);
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