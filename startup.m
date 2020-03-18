%%
%  file:   startup.m
%  author: Polcz Péter <ppolcz@gmail.com>
%
%  Created on Thu Aug 14 19:18:23 CEST 2014
%  Modified on 2019. June 05. (2019a)
%

global SCOPE_DEPTH
SCOPE_DEPTH = 0;

%% 

ROOT = cd;
EXTERNAL_DIR = [ROOT filesep 'results'];
EXTERNAL_LOG_DIR = [EXTERNAL_DIR filesep 'log'];

addpath(genpath([ROOT filesep 'ftools']))
addpath(genpath([ROOT filesep 'toolboxes']))
addpath(genpath([ROOT filesep 'utilities']))

EXTERNAL_DIR = cleanpath(EXTERNAL_DIR);
EXTERNAL_LOG_DIR = cleanpath(EXTERNAL_LOG_DIR);

setenv('ROOT',ROOT)
setenv('EXTERNAL_DIR',EXTERNAL_DIR);
setenv('EXTERNAL_LOG_DIR',EXTERNAL_LOG_DIR);

%% beginning of the scope

TMP_QVgVGfoCXYiYXzPhvVPX = pcz_dispFunctionName;

tbxmanager restorepath
disp 'Toolbox "SDPT3" in toolboxes/SDPT3-4.0 added to the Matlab path.'

P_init

pcz_dispFunctionEnd(TMP_QVgVGfoCXYiYXzPhvVPX);
clear TMP_QVgVGfoCXYiYXzPhvVPX
