function [ret] = snp_begin_end_scope(doch,event)
%%
%
%  file:   snp_begin_end_scope.m
%  author: Polcz Péter <ppolcz@gmail.com>
%
%  Created on 2016.01.31. Sunday, 15:04:58
%

[a,b] = pcz_generateBeginEndTimer;

try
    active = matlab.desktop.editor.getActive;
    editor = active.JavaEditor;
    
    linenr = editor.getLineNumber+1;
    lbegin = editor.lineAndColumnToPosition(linenr,0);
    editor.setCaretPosition(lbegin);
    
    b = strrep(b, '%%', '');
    b = regexprep(b, 'clear .*', '');
    b = strrep(b, newline, '');
    
    editor.insertTextAtCaret(a);
    clipboard('copy',b);
    
    % editor.appendText(sprintf('%s\n', b));
catch ex
    getReport(ex)
end


end
