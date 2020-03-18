function [ret] = snp_surround_with_key(doch,event)
%% 
%  
%  file:   snp_surround_with_key.m
%  author: Polcz Péter <ppolcz@gmail.com> 
%  
%  Created on 2016.02.01. Monday, 12:04:32
%

keyCode = event.getKeyCode()
shiftDown = event.isShiftDown()  % Inherited from java.awt.event.InputEvent

% Java objects concerning the active document and the editor
active = matlab.desktop.editor.getActive;
editor = active.JavaEditor;

% get selected text
selection = char(editor.getSelection);

cbegin = '';
cend = '';

switch keyCode
    case 91
        if shiftDown
            cbegin = '{';
            cend = '}';
        else
            cbegin = '[';
            cend = ']';
        end
            
end

editor.insertTextAtCaret([cbegin selection cend])

getID = event.getID()  % Inherited from java.awt.AWTEvent
getKeyChar = event.getKeyChar()
getKeyLocation = event.getKeyLocation()
getModifiers = event.getModifiers()  % Inherited from java.awt.event.InputEvent
getModifiersEx = event.getModifiersEx()  % Inherited from java.awt.event.InputEvent
isActionKey = event.isActionKey()
isAltDown = event.isAltDown()  % Inherited from java.awt.event.InputEvent
isAltGraphDown = event.isAltGraphDown()  % Inherited from java.awt.event.InputEvent
isConsumed = event.isConsumed()  % Inherited from java.awt.event.InputEvent
isControlDown = event.isControlDown()  % Inherited from java.awt.event.InputEvent
isMetaDown = event.isMetaDown()  % Inherited from java.awt.event.InputEvent
paramString = event.paramString()

% long getWhen()  % Inherited from java.awt.event.InputEvent

end