function [ret] = snp_copylinedown(doch,event)
%% 
%  
%  file:   snp_copylinedown.m
%  author: Polcz Péter <ppolcz@gmail.com> 
%  
%  Created on 2016.02.01. Monday, 12:03:06
%

%% duplicate line downwords

% Java objects concerning the active document and the editor
active = matlab.desktop.editor.getActive;
editor = active.JavaEditor;

selection = char(editor.getSelection);
cursor = editor.getCaretPosition;

if isempty(selection)
    linenr = editor.getLineNumber+1;

    lbegin = editor.lineAndColumnToPosition(linenr,0);
    lend = editor.lineAndColumnToPosition(linenr+1,0);
    lstr = editor.getText.substring(lbegin,lend);

    editor.setCaretPosition(lend)
    editor.insertTextAtCaret(lstr)
    editor.setCaretPosition(lend - lbegin + cursor)
else
    editor.insertTextAtCaret([selection selection]);
    editor.setSelection(editor.getCaretPosition - length(selection),editor.getCaretPosition);
end

end