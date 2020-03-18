function [return_val] = pcz_symzero(z, prec, N, varargin)
%% Script pcz_symzero
%  
%  file:   pcz_symzero.m
%  author: Peter Polcz <ppolcz@gmail.com> 
%  
%  Created on 2017.07.31. Monday, 23:54:28
%
%%

if nargin < 3 || isempty(N) || ischar(N)
    if nargin >= 3 && ischar(N)
        varargin = [N varargin];
    end
    N = 10;
end

if nargin < 2 || isempty(prec) || ischar(prec)
    if nargin >= 2 && ischar(prec)
        varargin = [prec varargin];
    end
    prec = 10;
end

ret = false;
if pcz_symzero1(z,prec,N,varargin{:})
    ret = true;
end

if ~ret
    s = symvar(sym(z));
    
    if ~isempty(s)    
        try
            [Theta,z0,q] = P_Pi_canonical_decomp(z(:), s);
            ret = pcz_symzero(Theta, prec, N);
        catch
            ret = false;
        end
    end
end

if nargout > 0
    return_val = ret;
else
    pcz_info(ret,varargin{:})
end

end