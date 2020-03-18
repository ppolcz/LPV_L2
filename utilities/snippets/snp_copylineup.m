function [ret] = snp_copylineup(doch,event)
%% 
%  
%  file:   snp_copylineup.m
%  author: Polcz Péter <ppolcz@gmail.com> 
%  
%  Created on 2016.02.01. Monday, 12:03:10
%

%% duplicate line upwords

% Java objects concerning the active document and the editor
active = matlab.desktop.editor.getActive;
editor = active.JavaEditor;

linenr = editor.getLineNumber+1;

lbegin = editor.lineAndColumnToPosition(linenr,0);
lend = editor.lineAndColumnToPosition(linenr+1,0);
cursor = editor.getCaretPosition;
lstr = editor.getText.substring(lbegin,lend);

editor.setCaretPosition(lbegin)
editor.insertTextAtCaret(lstr)
editor.setCaretPosition(cursor)

end