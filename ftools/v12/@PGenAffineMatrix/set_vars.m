function N = set_vars(N, vars_)
%% set_subsvars
%  
%  File: set_subsvars.m
%  Directory: 1_PhD_projects/00_my_toolboxes/algo_P/@PGenAffineMatrix
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2019.09.08. (szeptember  8, vasárnap), 13:07 [ATIRVA - OK]
%

%%

    vars_channels = symvar(N.channels);
    
    nem_szerepel = setdiff(vars_channels,vars_);
    
    if ~isempty(nem_szerepel)
        error('PGenAffineMatrix: Variables `%s` must be included!', ...
            strjoin(cellfun(@(x) {char(x)}, num2cell(nem_szerepel)),', '))
    end
    
    N.vars = vars_;
    N.xi = matlabFunction(N.channels,'vars',{N.vars});
    
    
end