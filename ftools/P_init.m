function [ret] = P_init(version)
%% 
%  
%  File: P_init.m
%  Directory: 1_PhD_projects/00_my_toolboxes/FinslerTools
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2018. September 26.
%  Modified on 2018. October 22.
%

%%

fname = mfilename('fullpath');
[path,~,~] = fileparts(fname);

% Remove every other version from path
v_all = dir([ path '/v*' ]);
for i = 1:numel(v_all)
    vi = v_all(i);
    gp = genpath([vi.folder '/' vi.name]);
    
    warning off MATLAB:rmpath:DirNotFound
    rmpath(gp)
    warning on MATLAB:rmpath:DirNotFound
end


% If no argument is given, choose the most recent version
if nargin < 1
    version = v_all(end).name;
end

% Numeric version argument to string
if isnumeric(version)
    if version < 10
        version = round(version * 10);
    end
    
    version = ['v' num2str(version)];
end


gp = genpath([path '/' version]);
addpath(gp)

end