function [str] = pcz_struct_append(str, varargin)
%% 
%  
%  file:   pcz_struct_append.m
%  author: Polcz Péter <ppolcz@gmail.com> 
%  
%  Created on 2017.01.24. Tuesday, 14:25:01
%
% Examples:
%  model = pcz_struct_append(model,vars,A,B,G,F,Pi, 'rewrite', true);
%  model = pcz_struct_append(model,[1 2 3],'name','v')
% 

assert(isstruct(str),'The first argument must be a struct');

props.snapshot = [ '_old_' pcz_fancyDate('var') ];
props.name = '';
props.rewrite = false;

%% SINGLE 
% This is for the case when we want to append a single, not certainly
% anonime variable
if      ... the third argument is a property name
        nargin > 2 && isempty(inputname(3)) && ischar(varargin{3}) ... 
        ... the number of arguments is odd (str, var, 'name','x','prop','value',...)
        && mod(nargin,2) == 0 && nargin >= 4
    
    props = parsepropval(props, varargin{2:end});
    
    if isvarname(props.name)
        append(props.name,varargin{1});
        return
    end
end

%% MULTIPLE
% This is for the case when we want to append several variables from
% existent in the workspace
opts_start = nargin;
for i = 2:nargin
    var = varargin{i-1};
    name = inputname(i);
    if isempty(name) && ischar(var) && ~isempty(var)
        opts_start = i - 1;
        break
    end
end

props = parsepropval(props, varargin{opts_start:end});

for i = 2:opts_start
    var = varargin{i-1};
    name = inputname(i);

    % fprintf('name = %s\n', name)
    % pcz_info(isempty(name), 'isempty(name)')
    % pcz_info(opts.rewrite, 'opts.rewrite')
    % pcz_info(isfield(str,name), 'isfield(str,name)')
    % pcz_info(~isempty(name) && (opts.rewrite || ~isfield(str,name)), '~isempty(name) && (opts.rewrite || ~isfield(str,name))')
    
    if ~isempty(name)
        append(name,var);
    end
end

    function append(name, value)
        if isfield(str, name) && ~check(name,value)
            newname = [ name props.snapshot ];
            str.(newname) = str.(name);
            warning('Field %s already exists, renamed to %s!', name, newname);
        end
        str.(name) = value;
    end

    function ret = check(name, value)
        oldval = str.(name);
        
        c(1) = all(size(oldval) == size(value));
        % c(2) = pcz_symeq(symvar(oldval),symvar(value));
        
        ret = all(c);
    end

end