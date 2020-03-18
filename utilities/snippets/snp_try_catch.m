function [ret] = snp_try_catch(doch,event)
%%
%
%  file:   snp_try_catch.m
%  author: Polcz Péter <ppolcz@gmail.com>
%
%  Created on 2016.01.31. Sunday, 18:46:04
%

newline = sprintf('\n');
tab = '    ';

try
    % Java objects concerning the active document and the editor
    active = matlab.desktop.editor.getActive;
    editor = active.JavaEditor;

    % get selected text
    selection = char(editor.getSelection);
    
    linenr = editor.getLineNumber;
    content = char(editor.getText);
    
    lines = strsplit(strrep(content,newline, [' ' newline]), newline);
    line = lines{linenr+1};
    
    % detect leading spaces in the actual line
    i = 1; while i <= numel(line) && line(i) == ' ', i = i+1; end
    lind = i-1;
    ind = repmat(' ', [1 lind]);
    
    % increase the indentation in the selected text after the linebreaks
    selection = strrep(selection, [newline ind], [newline ind tab]);
    
    % detect and trim the leading spaces in the selection
    i = 1; while i <= numel(selection) && selection(i) == ' ', i = i+1; end
    selection = selection(i:end);
    
    if isempty(selection)
        text = '';
    else
        text = sprintf([
            ind 'try\n' ...
            ind tab '%s\n' ...
            ind 'catch ex\n' ...
            ind tab 'disp(getReport(ex))\n' ...
            ind 'end\n' ...
            ], selection);
    end
    
    editor.insertTextAtCaret(text)

    %%
    
    linenr = editor.getLineNumber;
    lbegin = editor.lineAndColumnToPosition(linenr+1,0)
    lend = editor.lineAndColumnToPosition(linenr+2,0)
    editor.setSelection(lbegin,lend)
    editor.replaceText('', lbegin, lend)
%     editor.goToLine(linenr,0)
%     editor.insertTextAtCaret(num2str(l))
    newpos = editor.lineAndColumnToPosition(linenr+1,100000)
    editor.getCaretPosition
%     editor.setCaretPosition(newpos)
    
%     get(editor)
%     methods(editor, '-full')
%     
%     methods(editor.getDocument, '-full')
%     doc = editor.getDocument;
%     doctext = doc.getText
%     methods(doctext, '-full')
    
%     editor.replaceText('A',2,3)
    
catch ex
    getReport(ex)
end

ret = '';

end


% void appendText(java.lang.String)
% int getCaretPosition()
% int getLength()
% int getLineNumber()
% java.lang.String getLongName()
% java.lang.String getSelection()
% java.lang.String getShortName()
% java.lang.String getText()
% void goToLine(int,int)
% void goToLine(int,boolean)
% void goToPositionAndHighlight(int,int)
% void insertAndFormatTextAtCaret(java.lang.String)
% void insertTextAtCaret(java.lang.String)
% int lineAndColumnToPosition(int,int)
% int[] positionToLineAndColumn(int)
% void putProperty(java.lang.String,java.lang.Object)
% void refreshMenus()
% void reload()
% java.lang.String reloadAndReturnError()
% void removeEventListener(com.mathworks.matlab.api.editor.EditorEventListener)
% void removePropertyChangeListener(java.beans.PropertyChangeListener)
% void removePropertyChangeListener(java.lang.String,java.beans.PropertyChangeListener)
% void replaceText(java.lang.String,int,int)
% void save() throws java.lang.Exception
% java.lang.String saveAndReturnError()
% void saveAs(java.lang.String) throws java.lang.Exception
% java.lang.String saveAsAndReturnError(java.lang.String)
% void setCaretPosition(int)
% void setClean()
% void setDirtyUntilSave()
% void setEditable(boolean)
% void setEditorStatusBarText(java.lang.String)
% void setSelection(int,int)
% void setStatusText(java.lang.String)
% void smartIndentContents()
% java.lang.String toString()  % Inherited from java.lang.Object
% void unlock()
% void wait(long) throws java.lang.InterruptedException  % Inherited from java.lang.Object
% void wait(long,int) throws java.lang.InterruptedException  % Inherited from java.lang.Object
% void wait() throws java.lang.InterruptedException  % Inherited from java.lang.Object
