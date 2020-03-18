function [ret] = pcz_struct2nvpair(s)
%% pcz_struct2nvpair
%  
%  File: pcz_struct2nvpair.m
%  Directory: 2_demonstrations/lib/matlab
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2019. February 17.
%

%%

tmp = [fieldnames(s) struct2cell(s)].';
ret = tmp(:).';


end