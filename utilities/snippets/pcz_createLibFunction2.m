function [ret] = pcz_createLibFunction2(fn, varargin)
%% Script pcz_createLibFunction2
%  
%  File: pcz_createLibFunction2.m
%  Directory: include/matlab
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2018. February 01.
%

args = '';
if ~isempty(varargin)
    args = ['(' strjoin(varargin, ', ') ')'];
end

prefix = '';
if ~strcmp(fn(1:length(prefix)), prefix)
    fn = [prefix fn];
end

f = pcz_resolvePath([proot 'include/matlab/' fn]);

header = pcz_headerComment(f);
text = pcz_escape_sprintf(header);

fid = fopen(f.path, 'w');
fprintf(fid, ['function [ret] = ' f.bname args '\n' text '\n\nend']);

edit(f.path);

end