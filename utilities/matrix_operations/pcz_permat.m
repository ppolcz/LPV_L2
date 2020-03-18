function [I_sigma] = pcz_permat(sigma, sigma_good)
% I_sigma = pcz_permat(sigma)
% 
% Returns a permutation matrix I_sigma.
% 
% In each column and in each row only one nonzero element is present,
% which is 1. In the ith column the sigma(i)th element is 1.
% 
% Tensor representation: I_sigma(i,j) = delta(i,sigma(j))
% 
%  file:   pcz_permat.m
%  author: Peter Polcz <ppolcz@gmail.com> 
%  
%  Created on 2017.07.06. Thursday, 20:09:04
%

%%

N = numel(sigma);

I_sigma = eye(N);
if nargin >= 2
    I_sigma = I_sigma(sigma_good,:);    
end

I_sigma = I_sigma(:,sigma);


    
end