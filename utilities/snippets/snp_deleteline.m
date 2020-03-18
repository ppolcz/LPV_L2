function [ret] = snp_deleteline(doch,event)
%% 
%  
%  file:   snp_deleteline.m
%  author: Polcz Péter <ppolcz@gmail.com> 
%  
%  Created on 2016.02.01. Monday, 12:02:56
%

%% delete line

% Java objects concerning the active document and the editor
active = matlab.desktop.editor.getActive;
editor = active.JavaEditor;

linenr = editor.getLineNumber+1;
lbegin = editor.lineAndColumnToPosition(linenr,0);
lend = editor.lineAndColumnToPosition(linenr+1,0);

editor.replaceText('', lbegin, lend)


end