function pcz_struct_split(varargin)
%% 
%  
%  file:   pcz_struct_split.m
%  author: Polcz Péter <ppolcz@gmail.com> 
%  
%  Created on 2017.02.01. Wednesday, 20:46:41
%

first_prop = 1;
while nargin >= first_prop && isstruct(varargin{first_prop})
    first_prop = first_prop + 1;
end

props.snapshot = [ '_old_' pcz_fancyDate('var') ];
props.rewrite = true;
props = parsepropval(props, varargin{first_prop:end});

exclude = [ 
    pcz_var_exclude_patterns
    ];

exactly = @(name) ['^' name '$'];

for k = 1:numel(varargin)
    str = varargin{k};
    fn = fieldnames(str);

    if isvarname(inputname(k))
        exclude_all = [
            exclude
            exactly(inputname(k))
            ];
    end
    
    blacklist = pcz_regexp_match_bool(fn, exclude_all);

    % DEBUG = [fn num2cell(blacklist')]'

    for i = find(blacklist == 0)
        % disp([ fn{i} ' = '])
        % disp(str.(fn{i}))

        if evalin('caller', ['exist(''' fn{i} ''',''var'')']) ...
            && ~check(str.(fn{i}),evalin('caller',fn{i}))
            if props.rewrite
                newname = [ fn{i} props.snapshot ];
                assignin('caller', newname, evalin('caller',fn{i}))
                assignin('caller', fn{i}, str.(fn{i}))
                warning('Variable %s already exists, renamed to %s!', fn{i}, newname);
            else
                warning('Variable %s already exists, skipped!', fn{i}, newname);
            end
        else
            assignin('caller', fn{i}, str.(fn{i}))
        end
    end
end

    function ret = check(value, oldval)        
        c(1) = all(size(oldval) == size(value));
        % c(2) = pcz_symeq(symvar(oldval),symvar(value));
        
        ret = all(c);
    end

end