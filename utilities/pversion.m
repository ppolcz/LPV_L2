function [ret] = pversion
%% 
%  
%  file:   pversion.m
%  author: Polcz Péter <ppolcz@gmail.com> 
%  
%  Created on 2016.02.08. Monday, 15:40:56
%

release = version('-release');
ret = str2double(release(1:4)) + (release(5) == 'b')/2;


end