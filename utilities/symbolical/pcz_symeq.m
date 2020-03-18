function [ret] = pcz_symeq(a,b,varargin)
%% Script pcz_symeq
%  
%  file:   pcz_symeq.m
%  author: Peter Polcz <ppolcz@gmail.com> 
%  
%  Created on 2017.07.20. Thursday, 09:46:28
%
%%

z = a - b;
if isnumeric(z) || isempty(symvar(z))
    z = double(z);
else
    z = simplify(z);
end

if nargout > 0
    ret = pcz_symzero(z,varargin{:});
else
    pcz_symzero(z,varargin{:})
end

end