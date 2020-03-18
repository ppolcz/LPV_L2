function ret = pcz_generateSymStateVector(dim, name, assumpt, verb, be_null)
%%
%
%  file:   pcz_generateSymStateVector.m
%  author: Polcz Péter <ppolcz@gmail.com>
%
%  Created on 2016.01.09. Saturday, 19:45:02
%  Modified on 2019. June 10. (2019a)
%
%  Output: a command which generates x1, .. xdim, x = [x1, .. , xdim],
%  x_cell, x_str, x_vars
%
%%

% syms xi alpha real
% dim = 2;
% name = 'd';
% assumpt = [xi alpha];
% nargin = 3

if nargin < 1
    dim = 2;
end
nrs = num2cell(1:dim);


if nargin < 2
    name = 'x';
end

if nargin < 3
    assumpt = 'real';
end

if nargin < 4
    verb = 0;
end

if nargin < 5
    be_null = 0;
end


% { nargin, dim, name }

if nargin >= 3 ...
        && isa(assumpt,'sym') && isnumeric(dim) && isscalar(dim) ...
        && dim == numel(assumpt) && ischar(name)
    %% New interface
    symvars = cellfun(@(s) {char(s)}, num2cell(assumpt));
    dvars = cellfun(@(c)[name num2str(c)], nrs, 'uni', false);
    
    sym_stm = cellfun(@(var,def) {[var ' = ' def ';']}, dvars, symvars);
    sym_stm = [ strjoin(sym_stm, ' ') '\n' ];
else
    %% Old interface
    symvars = cellfun(@(c)[name num2str(c)], nrs, 'uni', false);
    dvars = symvars;
    
    if dim > 0
        sym_stm = ['syms ' strjoin(symvars, ' ') ' ' assumpt '\n'];
    elseif dim == 0
        sym_stm = '';
    end
end


%% Old core of the function
% Need to be defined: dvars, symvars, sym_stm

qsymvars = cellfun(@(c)['''' c ''''], symvars, 'uni', false);

x = name;
x_cell = [name '_cell'];
x_str = [name '_str'];
x_names = [name '_names'];
x_n = [name '_n'];
n_x = ['n_' name];
x_vars = [name '_vars'];

if ~be_null
    vars = {x x_cell x_str x_names x_n x_vars};
else
    vars = {x x_n x_vars};
end
qvars = cellfun(@(c)['''' c ''''], vars, 'uni', false);

if ~be_null
    text = [ sym_stm...
        x ' = [ ' strjoin(dvars, ' ; ') ' ];\n' ...
        x_cell ' = num2cell(' x ');\n' ...
        x_str ' = ''(' strjoin(symvars, ',') ')'';\n' ...
        x_names ' = { ' strjoin(qsymvars, ',') ' };\n' ...
        x_n ' = ' num2str(dim) ';\n' ...
        n_x ' = ' num2str(dim) ';\n' ...
        x_vars ' = { ' strjoin(qvars, ',') ' };\n' ];
else
    text = [ 
        x ' = zeros(' num2str(dim) ',1);\n' ...
        x_n ' = ' num2str(dim) ';\n' ...
        n_x ' = ' num2str(dim) ';\n' ...
        x_vars ' = { ' strjoin(qvars, ',') ' };\n' ];
end

command = sprintf(text);

if nargout > 0
    ret = command;
else
    if verb
        disp(command)
    end
    evalin('caller', command)
end

end
