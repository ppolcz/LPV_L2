function [ret] = snp_load_mat_file(doch,event)
%% 
%  
%  file:   snp_load_mat_file.m
%  author: Polcz Péter <ppolcz@gmail.com> 
%  
%  Created on 2016.06.07. Tuesday, 06:40:21
%

G_ = pglobals;

% Have user browse for a file, from a specified "starting folder."
% For convenience in browsing, set a starting folder from which to browse.
startingFolder = sprintf('%s/%s', proot, G_.RELPATH_GENERATED_FILES);
if ~exist(startingFolder, 'dir')
	% If that folder doesn't exist, just start in the current folder.
	startingFolder = pwd;
end
% Get the name of the mat file that the user wants to use.
defaultFileName = {
    fullfile(startingFolder, '*.mat;*.fig')
%     fullfile(startingFolder, '*.fig') 
    };
[baseFileName, folder] = uigetfile(defaultFileName, 'Select a mat file');
if baseFileName == 0
	% User clicked the Cancel button.
	return;
end
fullFileName = fullfile(folder, baseFileName);

[~,~,ext] = fileparts(fullFileName); 

if strcmp(ext,'.mat')
    load(fullFileName);
elseif strcmp(ext,'.fig');
    open(fullFileName)
end

end