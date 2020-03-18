function ret = pcs(varargin)
%% Script pcs
%  
%  file:   pcs.m
%  author: Peter Polcz <ppolcz@gmail.com> 
%  
%  Created on 2017. September 20.
%
%%

if nargout > 0
    ret = pcz_createScript(varargin{:});
else
    pcz_createScript(varargin{:});
end

end