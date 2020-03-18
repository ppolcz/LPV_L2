function cleaned_path_name = cleanpath(path_name, working_path)
%CLEANPATH Clean a file path name.
%
%   CLEANPATH(PATH_NAME, WORKING_PATH) removes redundant characters from
%   the path name, e.g. '//', '/./' as well as initial './' and circular
%   paths e.g. 'abc/def/../def/'.
%   Additionally, all file separators are set to the platform file
%   separator.
%
%   If a working path WORKING_PATH is given (the pwd function can be used
%   to get the current directory), PATH_NAME is assumed to be an absolute
%   path and WORKING_PATH is removed thus returning the relative path.
%

%
%   Marco Bucci rev. 2015/09/01
%	marco.bucci@infineon.com
%


nargs = nargin;

if (nargs < 2)
    working_path = [];
end;


% get the possible file separators
allowed_file_separators = '/';  % add the "standard" file separator
if filesep ~= '/' % add the system file separator if different
    allowed_file_separators = [filesep allowed_file_separators];
end


% remove insignificant (i.e. initial and ending) white spaces
path_name = strtrim(path_name);

% replace "wrong" file separators
for i=length(allowed_file_separators)
    path_name = strrep(path_name, allowed_file_separators(i), filesep);
end


% remove redundant './'
path_name = strrep(path_name, [filesep '.' filesep], filesep);

% remove double slashes
prev_len = 0;
while prev_len ~= length(path_name)
   prev_len = length(path_name);
   path_name = strrep(path_name, [filesep filesep], filesep);
end

% remove initial './' if path name is not empty
if length(path_name) > 2 && strcmp(path_name(1:2), ['.' filesep])
    path_name = path_name(3:end);
end

% remove redundant '../'
path_name = reduce_updir(path_name);

% get the relative path if requested
if ~isempty(working_path)
    path_name = relativepath(path_name, working_path);
end

cleaned_path_name = path_name;

end


function new_path_name = reduce_updir(path_name)
pre_updir ='';
[pre_updir, path_name] = r_reduce_updir(pre_updir, path_name);
new_path_name = [pre_updir path_name];
end

    
function [new_pre_updir, new_path_name] = r_reduce_updir(pre_updir, path_name)

while (length(path_name) > 3) && strcmp(path_name(1:3), ['..' filesep])
    pre_updir = [pre_updir '..' filesep];
    path_name = path_name(4:end);
end

% search for the first occurrence of '/../'
idx_2 = min(strfind(path_name, [filesep '..' filesep]));
if ~isempty(idx_2) % if found
    % search for the previous occurrence of '/'
    idx_1 = max(strfind(path_name(1:idx_2-1), filesep));
    path_name = [path_name(1:idx_1) path_name(idx_2+4:end)];
    % search for the next '../'to be removed
    [pre_updir, path_name] = r_reduce_updir(pre_updir, path_name);
end

new_pre_updir = pre_updir;
new_path_name = path_name;

end


function [path] = relativepath(absolute_path, working_path)

% clean the paths before comparing (call cleanpath itself, but without
% the working_path parameter)
absolute_path = cleanpath(absolute_path);
% save the last file separator if present
if absolute_path(end)==filesep
    last_sep = true;
else
    last_sep = false;
    % add a final file separator
    absolute_path(end+1) = filesep;
end

% clean the working directory path and, eventually, add a final file
% separator
working_path = cleanpath([working_path filesep]);

% search for a possible common path
min_len = min(length(absolute_path), length(working_path));
i = 0;
while i < min_len ...
        && absolute_path(i+1)==working_path(i+1)
    i = i+1;
end
last_idx = i;  % index of the end of the common path
% The common path must end with a file separator. Search for the last
% file separator inside the common path.
last_idx = max(strfind(working_path(1:last_idx), filesep));

% remove the common path
path = absolute_path(last_idx+1:end);
% remove the last file separator if it was not present
if ~last_sep
    path = path(1:end-1);
end

% search for the possible up directory '../' to be added
if last_idx ~= 0  % there is a common path
    s = '';
    % search for file separators in the not common part of the working
    % directory
    if last_idx < length(working_path)
        for i=length(working_path)-1:-1:last_idx
            if working_path(i) == filesep
                s = [s '..' filesep];
            end
        end
    end
    path = [s path];
end

if isempty(path)
    path = ['.' filesep];
end

end


