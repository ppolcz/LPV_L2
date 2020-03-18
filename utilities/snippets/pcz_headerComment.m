function [varargout] = pcz_headerComment(varargin)
%%
%
%  file:   pcz_headerComment.m
%  author: Polcz Peter (ppolcz@gmail.com)
%
%  Created on 2016.01.31. Sunday, 12:32:05
%  Reviewed on 2017. October 18. [to mlx]
% 

editor = 0;
fn = '';
bn = '';
dir = '';
nargin
if nargin > 0 && ischar(varargin{1})
    fn = varargin{1};
    [reldir,bn,~] = fileparts(fn);
    dir = strrep(pwd,proot,'');
    if ~isempty(reldir)
        dir = [ dir '/' reldir ];
    end
elseif nargin > 0 && isstruct(varargin{1})
    f = varargin{1};
    fn = f.fn;
    bn = f.bname;
    dir = f.reldir;
elseif nargin == 0
    active = matlab.desktop.editor.getActive;
    f = pcz_resolvePath(active.Filename);
    fn = f.fn;
    bn = f.bname;
    dir = f.reldir;
    editor = 1;
end


ret = sprintf([...
    '%%%% ' bn '\n'...
    '%%  \n'...
    '%%  File: ' fn '\n'...
    '%%  Directory: ' dir '\n'...
    '%%  Author: Peter Polcz (ppolcz@gmail.com) \n'...
    '%%  \n'...
    '%%  Created on ' pcz_fancyDate('informative') ' (' version('-release') ')\n'...
    '%%\n\n%%%%\n']);

if editor
    ret_text = sprintf([...
        ' File: ' fn '\n'...
        ' Directory: ' dir '\n'...
        ' Author: Peter Polcz (ppolcz@gmail.com) \n\n'...
        ' Created on ' pcz_fancyDate('informative') ' (' version('-release') ')\n']);
else
    ret_text = sprintf([...
        'File: ' fn '\n'...
        'Directory: ' dir '\n'...
        'Author: Peter Polcz (ppolcz@gmail.com) \n\n'...
        'Created on ' pcz_fancyDate('informative') ' (' version('-release') ')\n']);
end
    
if nargout > 0
    varargout{1} = ret;
else
    clipboard('copy', ret_text)
end

end