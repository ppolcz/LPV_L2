function [ret] = Lie(f, g, x)
%% Lie
%  
%  File: Lie.m
%  Directory: 7_ftools/utilities/Lie
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2019. October 30. (2019a)
%

%%

narginchk(3,3)

namef = 'f';
if ~isempty(inputname(1))
    namef = inputname(1);
end

namex = 'x';
if ~isempty(inputname(3))
    namex = inputname(3);
end

if numel(x) ~= size(f,1)
    error('dim(%s) = %d, while size(%s,1) = %d', namex, numel(x), namef, size(f,1));
end

ret = jacobian(g,x) * f;

end