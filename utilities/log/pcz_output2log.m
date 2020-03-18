function [ret] = pcz_output2log(logname,varargin)
%% pcz_output2log
%
%  File: pcz_output2log.m
%  Directory: 2_demonstrations/lib/matlab
%  Author: Peter Polcz (ppolcz@gmail.com)
%
%  Created on 2019. February 28.
%

args.rewrite = 1;
args = parsepropval(args,varargin{:});

%%

% logname = '/home/ppolcz/_/users/public_html/files/results/CDC-2019-Ex1/2019-02-28_17-12-39_i1441_LPV_passivity_CDC2019_4D_3x2_p2.txt';

sed_commands = {
    's/<\/strong>|<strong>//g'
    's/<a href="[^<>]*">([^<>]*)<\/a>/\1/g'
    's/\].$//g'
    's/\[.INFO\]./INFO/g'
    's/(\[|\]|\}|\{)[^0-9a-zA-Z ]//g'
    's/Saved to workspace: .*//g'
    's/opentoline\(.([^,]*).,([0-9]*),.*\)/\1:\2/g'
    % 's/┌/|-/g'
    % 's/│/|/g'
    % 's/└/|-/g'
    };

[path,name,ext] = fileparts(logname);

new_logname = sprintf('%s/%s_output%s',path,name,ext);

sed_commands = strjoin(sed_commands,';');

cmd = sprintf('cat %s | sed -r ''%s'' > %s', logname, sed_commands, new_logname);

output = system(cmd);

if args.rewrite == 1
    movefile(new_logname,logname);    
end


% pcz_dispFunction(output)

end