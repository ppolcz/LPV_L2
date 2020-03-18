function [ret] = pcz_sym_full(name,dim1,dim2,startindex)
%% pcz_sym_full
%  
%  File: pcz_sym_full.m
%  Directory: 2_demonstrations/lib/matlab
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2018. October 29.
%

%%

if nargin < 4
    startindex = 1;
end

syms = pcz_sym_indexed(name,dim1*dim2,startindex);

ret = reshape(syms, [dim1,dim2]);

end