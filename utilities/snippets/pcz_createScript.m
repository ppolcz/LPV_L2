function ret = pcz_createScript(fn)
%% 
%  File: pcz_createScript.m
%  Directory: 2_demonstrations/lib/matlab
%  Author: Peter Polcz (ppolcz@gmail.com) 
% 
%  Created on Thu Aug 14 19:03:34 CEST 2014
%  Modified on 2018. April 16.

%%

f = pcz_resolvePath(fn);
clear fn

[~,globscode] = pglobals;

globalScope = [
    sprintf('%% Automatically generated stuff\n') ...
    newline ...
    globscode ...
    ];

header = pcz_headerComment(f);

beginning = [
    newline ...
    sprintf('\n%s\n\n', snp_persistent) ...
    '%%' newline
    ];

text = [ header globalScope beginning 'persist.stoplog;' newline ];
text = pcz_escape_sprintf(text);

if nargout > 0
    ret = text;
    return
end

fid = fopen(f.path, 'w');
fprintf(fid, text);
edit(f.path);

end
