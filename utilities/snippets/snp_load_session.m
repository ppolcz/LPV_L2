function [ret] = snp_load_session(doch,event)
%% 
%  
%  file:   snp_load_session.m
%  author: Polcz Péter <ppolcz@gmail.com> 
%  
%  Created on 2016.06.29. Wednesday, 11:23:48
%


G_ = pglobals;

% Have user browse for a file, from a specified "starting folder."
% For convenience in browsing, set a starting folder from which to browse.
startingFolder = sprintf('%s/%s', proot, G_.RELPATH_SESSION);
if ~exist(startingFolder, 'dir')
	% If that folder doesn't exist, just start in the current folder.
	startingFolder = pwd;
end
% Get the name of the mat file that the user wants to use.
defaultFileName = {
    fullfile(startingFolder, '*.txt')
%     fullfile(startingFolder, '*.fig') 
    };
[baseFileName, folder] = uigetfile(defaultFileName, 'Select a mat file');
if baseFileName == 0
	% User clicked the Cancel button.
	return;
end
fullFileName = fullfile(folder, baseFileName);

fid = fopen(fullFileName, 'r');
tline = fgets(fid);
while ischar(tline)
    tline = strtrim(tline);
    fprintf('Opening |%s| ...\n', tline)
    if exist(tline,'file')
        open(tline)
    else
        warning('File `%s` not found', tline);
    end
    tline = fgets(fid);
end
fclose(fid);

end