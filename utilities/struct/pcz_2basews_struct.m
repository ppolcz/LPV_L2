function pcz_2basews_struct(varargin)
%% pcz_2basews_struct
%  
%  File: pcz_2basews_struct.m
%  Directory: utilities/struct
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2020. March 19. (2019b)
%

%%

timp = now; 
varname = sprintf('r%sh_%sm_%ss', ...
    datestr(timp, 'yyyy_mmmmdd_HH'), ...
    datestr(timp,'MM'), ...
    datestr(timp,'SS'));



if nargin == 0
    vars = evalin('caller', 'who');
    vars = vars(~strncmp(vars,'persist',7));

    n = numel(vars);
else
    vars = varargin;
    n = nargin;
end

s = struct;
for i = 1:n
    if i <= nargin && ~isempty(inputname(i))
        s.(inputname(i)) = vars{i};
        vars{i} = inputname(i);
    elseif isvarname(vars{i})
        s.(vars{i}) = evalin('caller', vars{i});
    end
end

RESULTS__.(varname) = s;
assignin('base','RESULTS__',RESULTS__);

vars = join(vars, ', ');

pcz_dispFunction('Saved to workspace.RESULTS__.%s: %s', varname, vars{1});

end