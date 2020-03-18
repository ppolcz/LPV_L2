function [ret] = pround(x, prec)
%% 
%  
%  file:   pround.m
%  author: Polcz Péter <ppolcz@gmail.com> 
%  
%  Created on 2016.04.19. Tuesday, 14:51:09
%

error "Use round instead"

multiplier = 10^prec;

ret = round(multiplier * x) / multiplier;

end