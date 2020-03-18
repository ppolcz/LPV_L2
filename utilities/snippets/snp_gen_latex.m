function [ret] = snp_gen_latex(doch,event)
%% 
%  
%  file:   snp_gen_latex.m
%  author: Polcz Péter <ppolcz@gmail.com> 
%  
%  Created on 2016.02.21. Sunday, 16:13:45
%
% @syms a b c d real
% @sed s/Kutya/K_{utya}/g
% 
% @syms a b c d e x1 real
% @align 2
% h = a+b + exp(b) + 10^b * b^2 * x1^2;
% g = e;
% k = [a x1 x1^2 ; x1 b c]; % asda
% l = [a x1 x1^2 ; x1 b c]; 
% o = [a [x1 x1^2] ; x1 b c];

title = 'Generate LaTeX Snippet';
separator = repmat('-', [1,length(title) + 6]);
newline = sprintf('\n');

pprintc('\n{1}\n-- {2} --\n{1}\n\nVariables defined:\n\n',separator, title)

% Java objects concerning the active document and the editor
active = matlab.desktop.editor.getActive;
editor = active.JavaEditor;
text = char(editor.getText);

sed_cmd = '';

syms mu_max real

m_syms = sprintf(['%%\\s*@(?<first>syms [^\n]*)|\n(?<first>syms [^\n]*)|' ... 
    '\n(?<first>eval.pcz_generateSymStateVector[^\n]*)']);
m_sed = '%\s*@sed';

matches = regexp(text, m_syms, 'names');
for i = 1:numel(matches)
    fprintf('%s\n', matches(i).first)
    eval(matches(i).first)
end

% TODO it would be more elegant if using regexp named tokens
sed_positions = regexp(text, m_sed);
for pos = sed_positions
    [linenr,~] = pos2linecol(editor, pos);
    line = getLine(editor, linenr);
    
    [from, to] = regexp(line, 's\/.*\/g');
    sed_cmd = [sed_cmd line(from:to) ';']; %#ok<AGROW>
end

no_brackets = '[^\[\]]';
simple_brackets = [ '\[' no_brackets '*\]' ];
double_brackets = [ '\[(' no_brackets '*(' simple_brackets ')*)*\]' ];
patterns = {
    % '% @align 2'
    '%\s*@align\s*(?<align>\d+)'

    % h = a+b + exp(b) + 10^b * b^2 * x1^2;
    '(?<varname>\w+)\s*=\s*(?<def>[^\[\];,]*)\s*[;,]'

    % l = [a x1 x1^2 ; x1 b c];
    [ '(?<varname>\w+)\s*=\s*(?<def>' simple_brackets ')[\s,;]*' ]
        
    % o = [a [x1 x1^2] ; x1 b c];
    [ '(?<varname>\w+)\s*=\s*(?<def>' double_brackets ')\s*' ]  
    }

selection = char(editor.getSelection());
matches = regexp(selection, strjoin(patterns', '|'), 'names');

pprintc('\n{1}\nLaTeX code:\n\n', separator);

align = 0;
for i = 1:numel(matches)
    
    if ~isempty(matches(i).align)
        align = str2double(matches(i).align);
        continue;
    end
    
    varname = matches(i).varname;
    try
        var = eval(matches(i).def);

        stage1 = pcz_latex_v16k(var, 'disp_mode', align, 'varname', varname);
        
        system_command = ['echo ''' stage1 ''' | sed -r "' sed_cmd '"'];
        [status,cmdout] = system(system_command);
        
        if ~status
            % fprintf('%% sed: %s\n%% %s\n', sed_cmd, system_command)
            fprintf('%s', cmdout)
        else
            fprintf('sed failed, error id: %d\nmessage = %s\ncommand: `%s`\n', status, cmdout, system_command);
            fprintf('%s', stage1)
        end
    catch ex
        ex.getReport
    end
end

disp ' '

function [linenr,colnr] = pos2linecol(editor, pos)
linecol = editor.positionToLineAndColumn(pos);
linenr = linecol(1);
colnr = linecol(2)-1;

function line = getLine(editor, linenr)
lbegin = editor.lineAndColumnToPosition(linenr,0);
lend = editor.lineAndColumnToPosition(linenr+1,0);
line = char(editor.getText.substring(lbegin,lend-1));

