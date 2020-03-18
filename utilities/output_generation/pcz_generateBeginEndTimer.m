function varargout = pcz_generateBeginEndTimer
%% 
%  
%  File: pcz_generateBeginEndTimer.m
%  Directory: 2_demonstrations/lib/matlab
%  Author: Peter Polcz (ppolcz@gmail.com) 
% 
%  Created on 2016.01.17. Sunday, 13:30:51
%  Modified on 2018. April 16.
%

%
%%

var = ['TMP_' pcz_generateString(20, 0) ];
beginning = sprintf('%s = pcz_dispFunctionName;\n', var);
ending = sprintf('\npcz_dispFunctionEnd(%s);\nclear TMP_*\n', var);

%%
if nargout == 1
    varargout{1} = [beginning ending];
elseif nargout == 2
    varargout{1} = beginning;
    varargout{2} = ending;
else
    beginning = sprintf('pcz_dispFunction('''')\n%s = pcz_dispFunctionName;\n', var);
    ending = sprintf('pcz_dispFunctionEnd(%s);', var);
    clipboard('copy', [beginning ending])
end
    
end
