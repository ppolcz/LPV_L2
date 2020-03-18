function [varargout] = pcz_createClass(fn)
%% 
%  
%  file:   pcz_createClass.m
%  author: Polcz Péter <ppolcz@gmail.com> 
% 
%  Created on Sat Aug 23 17:35:56 CEST 2014
%

%%

f = pcz_resolvePath(fn);

header = pcz_headerComment(f);

if exist(f.path, 'file') 
    edit(f.path)
    return
end

% MyDate = java.util.Date;
% StringDate = javaMethod('toString',MyDate);

fid = fopen(f.path, 'w');

fprintf(fid, ['classdef ' f.bname '\n']);
fprintf(fid, pcz_escape_sprintf(header));

fprintf(fid,'properties (Constant = true)\nend\n\n');
fprintf(fid,'properties (GetAccess = private, SetAccess = private)\nend\n\n');
fprintf(fid,'properties (GetAccess = public, SetAccess = public)\nend\n\n');
fprintf(fid,'methods (Access = public)\n%%%% public\nend\n\n');
fprintf(fid,'methods (Access = private)\n%%%% private\nend\n\n');
fprintf(fid, 'end\n');
        
edit(f.path);

end
