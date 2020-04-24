function [ret] = P_cartprod(A,B)
%% 
%  
%  file:   P_cartprod.m
%  author: Polcz Péter <ppolcz@gmail.com> 
%  
%  Created on 2016.04.19. Tuesday, 09:21:07
%

if isempty(B), B = zeros(1,0); end
if isempty(A), A = zeros(1,0); end

ret = [
    reshape(repmat(A', [size(B,1) 1]), [size(A,2), size(A,1)*size(B,1)])' ...
    repmat(B, [size(A,1), 1])
    ];

end