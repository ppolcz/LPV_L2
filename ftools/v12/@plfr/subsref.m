function varargout = subsref(pLFR, S)
%% subsref
%  
%  File: subsref.m
%  Directory: 1_PhD_projects/00_my_toolboxes/FinslerTools/v11/@plfr
%  Author: Peter Polcz (ppolcz@gmail.com) 
% 
%  Created on 2019. March 09.
% 

%%


switch S(1).type
    case '()'
        
        if numel(S.subs) == 1 && iscell(S.subs{1})
            x = S.subs;
        else
            % ([1;2],[3 4],5) --> {1 ; 2 ; 3 ; 4 ; 5}
            subs = cellfun(@(c) {c(:)},S.subs);
            x = num2cell(vertcat(subs{:}));
        end
                
        if numel(x) == numel(pLFR.subsvars)    
            
            % Regi megoldas
            % varargout{1} = pLFR.fh(x{:});
            
            % Uj megoldas: 2019.09.14. (szeptember 14, szombat), 19:06 
            % Egy nullaval kiegeszitve: 2019.11.06. (november  6, szerda), 01:34            
            d = pLFR.delta_fh(x{:},zeros(size(x{1})));

            varargout{1} = pLFR.A + pLFR.B/( pLFR.I - diag(d)*pLFR.D)*diag(d)*pLFR.C;

        elseif nargout > 0
            varargout = cell(1,nargout);
            [varargout{:}] = builtin('subsref', pLFR, S);
        else
            builtin('subsref', pLFR, S);
        end
                
    case '{}'
        error('{} indexing not supported');
    case '.'
        % isprop(N,S(1).subs)
        if any(strcmp(fieldnames(pLFR),S(1).subs))
            varargout{1} = builtin('subsref', pLFR, S);
        elseif nargout > 0
            varargout = cell(1,nargout);
            [varargout{:}] = builtin('subsref', pLFR, S);
        else
            builtin('subsref', pLFR, S);
        end
end

end