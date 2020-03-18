function [ret] = snp_pevalrep(doch,event)
%% 
%  
%  file:   snp_pevalrep.m
%  author: Polcz Péter <ppolcz@gmail.com> 
%  
%  Created on 2016.02.14. Sunday, 16:45:06
%

% Java objects concerning the active document and the editor
active = matlab.desktop.editor.getActive;
editor = active.JavaEditor;

% get selected text
selection = char(editor.getSelection);

if length(selection) < 6 || ~strcmp(selection(1:6),'peval ')
    selection = [ 'peval ' selection ];
end
selection = strrep([ selection ' #r ;' ], sprintf('\n'), ' ');

eval(selection); ret = ans;
editor.insertTextAtCaret(ret)


%% Test strings
% struct _ret kutyagumi alma korte '' ?farok torony #tex

end