function [str] = pcz_struct_append_wspvars(str,varargin)
%% 
%  
%  file:   pcz_struct_append_wspvars.m
%  author: Polcz Péter <ppolcz@gmail.com> 
%  
%  Created on 2017.02.01. Wednesday, 20:14:21
%
% Examples:
%  str = pcz_struct_append_wspvars;
%  str = pcz_struct_append_wspvars(str);
%  str = pcz_struct_append_wspvars('rewrite',false);

if nargin == 0 || ~isstruct(str)
    str = struct;
    strname = '[noname]';
    
    if nargin > 0 ; varargin = [str varargin]; end
else
    strname = inputname(1);
end

props.exclude = {};
props = parsepropval(props, varargin{:});

exactly = @(name) ['^' name '$'];
exclude = [ 
    pcz_var_exclude_patterns
    exactly(strname)
    props.exclude
    ];


vars = evalin('caller', 'who');

blacklist = pcz_regexp_match_bool(vars, exclude);

% DEBUG = [vars num2cell(blacklist')]'


for j = find(blacklist == 0)
    var = vars{j};
    % disp(var)
    % disp(evalin('caller',var))
    
    % do not append itself
    str = pcz_struct_append(str, evalin('caller',var), 'name', var, varargin{:});
end




end