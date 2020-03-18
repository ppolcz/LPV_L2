function [ret] = Lbr(f, g, x)
%% Lbr
%  
%  File: Lbr.m
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

nameg = 'g';
if ~isempty(inputname(2))
    nameg = inputname(2);
end

namex = 'x';
if ~isempty(inputname(3))
    namex = inputname(3);
end

if numel(x) ~= size(f,1)
    error('dim(%s) = %d, while size(%s,1) = %d', namex, numel(x), namef, size(f,1));
end

if numel(x) ~= size(g,1)
    error('dim(%s) = %d, while size(%s,1) = %d', namex, numel(x), nameg, size(g,1));
end

ret = Lie(f,g,x) - Lie(g,f,x);

end