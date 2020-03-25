function [ret] = pcz_dispFunction_scalar(varargin)
%% pcz_dispFunction_scalar
%  
%  File: pcz_dispFunction_scalar.m
%  Directory: utilities/output_generation
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2020. March 19. (2019b)
%

for i = 1:nargin
    pcz_dispFunction2('%s = %s', inputname(i), num2str(varargin{i}))
end

end