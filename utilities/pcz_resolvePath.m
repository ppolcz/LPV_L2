function ret = pcz_resolvePath(varargin)
%% 
%  
%  file:   pcz_resolvePath.m
%  author: Polcz Péter <ppolcz@gmail.com> 
%  
%  Created on 2016.01.31. Sunday, 12:40:36
%  Reviewed on 2016.02.19. Friday, 00:43:33
%
%

% initialize struct ret
ret.fn = '';
ret.bname = '';
ret.ext = '';
ret.reldir = '';
ret.relpath = '';
ret.dir = '';
ret.path = '';
ret.type = '';

ROOT = getenv('ROOT');

fn = varargin{1};
[reldir,ret.bname,ret.ext] = fileparts(fn);

if isempty(ret.bname), return, end
if isempty(ret.ext), ret.ext = '.m'; end

ret.fn = [ret.bname ret.ext];

if isempty(reldir) || exist([cd '/' reldir], 'dir')
    % relative directory to the current dir
    if ~isempty(reldir) && reldir(1) ~= '/', reldir = ['/' reldir]; end
    ret.type = 'cd';
    ret.dir = [cd reldir];
elseif exist([ROOT '/' reldir], 'dir')
    % relative directory to the workspace dir
    ret.type = 'wroot';
    ret.dir = [ROOT '/' reldir];
elseif exist(reldir,'dir')
    % relative directory to system root (i.e. absolute path)
    ret.type = 'root';
    ret.dir = reldir;
else
    % not found

    mkdir(reldir);
    ret.type = 'newly_created';
    ret.dir = reldir;

    % warning(['directory not found: ' reldir])
end

ret.path = [ret.dir '/' ret.fn];

ret.dir = char(java.io.File(ret.dir).getCanonicalPath());
ret.path = char(java.io.File(ret.path).getCanonicalPath());

ret.reldir = strrep(ret.dir, [ ROOT '/' ], '');
ret.relpath = strrep(ret.path, [ ROOT '/' ], '');

% ennek nem itt lenne a helye:
% if exist(ret.path, 'file')
%     edit(ret.path);
%     error(['File ' ret.path ' already exists!']);
% end

end