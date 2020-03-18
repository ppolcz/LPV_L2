function snp_function_header(doch,event)
%%
%
%  file:   snp_function_header.m
%  author: Polcz Péter <ppolcz@gmail.com>
%
%  Created on 2016.01.31. Sunday, 17:38:20
%

try
    active = matlab.desktop.editor.getActive;
    editor = active.JavaEditor;

    f = pcz_resolvePath(active.Filename);
    
    ret = sprintf('function [varargout] = %s(varargin)\n', f.bname);
    
    editor.setCaretPosition(0);
    editor.insertTextAtCaret(ret);
    editor.appendText(sprintf('\nend'));    
catch ex
    getReport(ex)
end

end

