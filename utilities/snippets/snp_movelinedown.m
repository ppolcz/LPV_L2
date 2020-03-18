function [ret] = snp_movelinedown(doch,event)
%% 
%  
%  file:   snp_movelinedown.m
%  author: Polcz Péter <ppolcz@gmail.com> 
%  
%  Created on 2016.02.01. Monday, 12:02:50
%

%% move line downwords

% Java objects concerning the active document and the editor
active = matlab.desktop.editor.getActive;
editor = active.JavaEditor;

linenr = editor.getLineNumber+1;

lbegin = editor.lineAndColumnToPosition(linenr,0);
lend = editor.lineAndColumnToPosition(linenr+1,0);
lend2 = editor.lineAndColumnToPosition(linenr+2,0);

if lend < lend2
    cursor = editor.getCaretPosition;
    lstr = char(editor.getText.substring(lbegin,lend));
    lstr2 = char(editor.getText.substring(lend,lend2));

    editor.replaceText([lstr2 lstr], lbegin, lend2);
    editor.setCaretPosition(lend2-lend+cursor)
end

end