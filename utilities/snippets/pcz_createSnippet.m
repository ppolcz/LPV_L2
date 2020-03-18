function [ret] = pcz_createSnippet(fn)
%% 
%  
%  file:   pcz_createSnippet.m
%  author: Polcz Péter <ppolcz@gmail.com> 
%  
%  Created on 2016.01.31. Sunday, 12:31:13
%

prefix = 'snp_';
if ~strcmp(fn(1:length(prefix)), prefix)
    fn = [prefix fn];
end

f = pcz_resolvePath(['include/matlab/snippets/' fn]);

header = pcz_headerComment(f);
text = pcz_escape_sprintf(header);

fid = fopen(f.path, 'w');
fprintf(fid, ['function [ret] = ' f.bname '(doch,event)\n' text '\n\nend']);

edit(f.path);

end