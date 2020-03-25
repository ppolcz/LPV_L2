function [ret] = pcz_2basews(varargin)
%% pcz_2basews
%  
%  File: pcz_2basews.m
%  Directory: 2_demonstrations/lib/matlab
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2018. November 21.
%

%%

if nargin == 0
    vars = evalin('caller', 'who');
    vars = vars(~strncmp(vars,'persist',7));

    n = numel(vars);
else
    vars = varargin;
    n = nargin;
end

for i = 1:n
    if i <= nargin && ~isempty(inputname(i))
        assignin('base',inputname(i),vars{i});
        vars{i} = inputname(i);
    elseif isvarname(vars{i})
        assignin('base',vars{i},evalin('caller', vars{i}));
    end
end

vars = join(vars, ', ');

pcz_dispFunction('Saved to workspace: %s', vars{1});

end