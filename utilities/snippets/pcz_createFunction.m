function pcz_createFunction(fn, varargin)
%% 
%  
%  file:   pcz_createFunction.m
%  author: Polcz Péter <ppolcz@gmail.com> 
% 
%  Created on Thu Aug 14 19:03:34 CEST 2014
%

%%

args = '';
if ~isempty(varargin)
    args = ['(' strjoin(varargin, ', ') ')'];
end

f = pcz_resolvePath(fn);

header = pcz_headerComment(f);
text = pcz_escape_sprintf(header);

fid = fopen(f.path, 'w');
fprintf(fid, ['function [ret] = ' f.bname args '\n' text '\n\nend']);

edit(f.path);

end