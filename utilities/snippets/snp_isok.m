function [ret] = snp_isok(doch,event)
%% 
%  
%  file:   snp_isok.m
%  author: Polcz Péter <ppolcz@gmail.com> 
%  
%  Created on 2016.02.14. Sunday, 16:44:50
%

shiftDown = event.isShiftDown();  % Inherited from java.awt.event.InputEvent

% Java objects concerning the active document and the editor
active = matlab.desktop.editor.getActive;
editor = active.JavaEditor;  

linenr = editor.getLineNumber+1;
lbegin = editor.lineAndColumnToPosition(linenr,0);
lend = editor.lineAndColumnToPosition(linenr+1,0)-1;

cursor = editor.getCaretPosition;

if shiftDown
    lstr = char(editor.getText.substring(lbegin,lend));
    lstr = regexprep(lstr, '%.*$', '');
    editor.replaceText(lstr, lbegin, lend);    
else
    editor.setCaretPosition(lend)
    editor.insertTextAtCaret(' %#ok')
    editor.setCaretPosition(cursor)
end
    
editor.setCaretPosition(cursor)
    
end