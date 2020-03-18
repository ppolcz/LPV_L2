function [ret] = pcz_sym_indexed(name, nr, startindex)
%% 
%  
%  File: pcz_sym_indexed.m
%  Directory: 2_demonstrations/lib/matlab
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2018. October 23.
%

%%

narginchk(2,3)

if nargin <= 2
    startindex = 1;
end

%%

assert(startindex >= 0, 'start index should be non-negative')

indices = startindex:nr+startindex-1;

assert(numel(indices) == nr);

minwidth = min(floor(log10(indices)) + 1);
maxwidth = max(floor(log10(indices)) + 1);

ret = sym(zeros(1,nr));

for i = minwidth:maxwidth
    tmp = sym([ name repmat('0',[1, maxwidth-i]) ],[10^i-1 , 1]);
    ami_kell = max(startindex,10^(i-1)):min(nr+startindex,10^i)-1;
    
    ret(ami_kell-startindex+1) = tmp(ami_kell);
end

end

%%

function self_check
%%

b = pcz_indexed_sym('b',8,6)
assert(all(size(b) == [1,8]))

a = pcz_indexed_sym('a',100,8);
assert(all(size(a) == [1,100]))
[a(1:2), a(end)]

a = pcz_indexed_sym('a',100,1001);
assert(all(size(a) == [1,100]))
[a(1:2), a(end)]

a = pcz_indexed_sym('a',100,989);
assert(all(size(a) == [1,100]))
[a(1:2), a(end)]

end