function snp_run_line
%% 
%  
%  file:   snp_run_line.m
%  author: Polcz Péter <ppolcz@gmail.com> 
%  
%  Created on 2016.11.20. Sunday, 15:56:36
%
% <https://www.mathworks.com/matlabcentral/answers/132119-keyboard-shortcut-to-evaluate-current-line
% Evluate current line>

currentEditor = matlab.desktop.editor.getActive;
originalSelection = currentEditor.Selection;

% Select the whole line
currentEditor.Selection = [originalSelection(1) 1 originalSelection(3) Inf];
selected = currentEditor.SelectedText;

while ~isempty(strfind(currentEditor.SelectedText,sprintf('...')))
    currentEditor.Selection = [currentEditor.Selection(3)+1 1 currentEditor.Selection(3)+1 Inf];
    selected = [selected, currentEditor.SelectedText];
end

% strrep(selected,'...','')

% remove ellipse and evaluate
evalin('caller', strrep(selected,'...','')); 

% Advance 1 line
currentEditor.Selection = [currentEditor.Selection(3)+1 1 currentEditor.Selection(3)+1 1]; 
clear currentEditor originalselection

end

function tmp
%%
fprintf('%s\n\n',...
    'kutya')

fprintf('%s\n\n','farka')

disp 'kutya ... asda'

end