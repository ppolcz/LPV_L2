function pcz_createLibFunction(fn, varargin)
%% 
%  
%  file:   pcz_createLibFunction.m
%  author: Polcz Péter <ppolcz@gmail.com> 
%  
%  Created on 2017.01.06. Friday, 12:53:50
%

args = '';
if ~isempty(varargin)
    args = ['(' strjoin(varargin, ', ') ')'];
end

prefix = '';
if ~strcmp(fn(1:length(prefix)), prefix)
    fn = [prefix fn];
end

f = pcz_resolvePath(['2_demonstrations/lib/matlab/' fn]);

header = pcz_headerComment(f);
text = pcz_escape_sprintf(header);

fid = fopen(f.path, 'w');
fprintf(fid, ['function [ret] = ' f.bname args '\n' text '\n\nend']);

edit(f.path);

end