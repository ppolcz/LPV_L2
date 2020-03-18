function [ret] = snp_varargin
%% 
%  
%  file:   snp_varargin.m
%  author: Polcz Péter <ppolcz@gmail.com> 
%  
%  Created on 2016.01.31. Sunday, 12:16:17
%
%% beginning of the scope
TMP_IbSWJNMuIiKbocfQKqXb = pcz_dispFunctionName;

%%

% fname = full path of the actual file
eval(pcz_cmd_fname('fname'));

% global variables
G__ = pglobals; 

% cache directory for the mat files
%persist_bin_dir = proot(fname,G__.RELPATH_MAT_FILES);


%% end of the scope
pcz_dispFunctionEnd(TMP_IbSWJNMuIiKbocfQKqXb);
clear TMP_IbSWJNMuIiKbocfQKqXb

end