%%
%  file:   startup.m
%  author: Polcz Péter <ppolcz@gmail.com>
%
%  Created on Thu Aug 14 19:18:23 CEST 2014
%  Modified on 2019. June 05. (2019a)
%  Major review on 2020. April 22. (2019b)
%

%% 

ROOT = cd;
A_MAT_PATH = [ROOT filesep 'A.mat'];

setenv('ROOT',ROOT)
setenv('A_MAT_PATH',A_MAT_PATH);

addpath(genpath('ftools'))
addpath(genpath('utilities'))

% cd ../../7_ftools
% tbxmanager restorepath
% addpath(genpath('toolboxes_other'))
% cd(ROOT)
