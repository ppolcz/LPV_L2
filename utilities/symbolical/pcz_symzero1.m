function [ret] = pcz_symzero1(z, prec, N, varargin)
%% Script pcz_symzero1
%  
%  file:   pcz_symzero1.m
%  author: Peter Polcz <ppolcz@gmail.com> 
%  
%  Created on 2017.07.31. Monday, 23:53:35
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

s = symvar(sym(z));

% prec
% N
% varargin

if isempty(s)
    ZERO = double(z);
    ret = all(abs(ZERO(:)) < 10^(-prec));
else
    z_fh = matlabFunction(z(:), 'vars', {s(:)});
    ZERO = zeros(numel(z),N);    
    for i=1:N
        ZERO(:,i) = z_fh(rand(numel(s),1));
    end
    
    greater = sum(abs(ZERO) > 10^(-prec),2);
    indices = find(greater);
    perc = numel(indices) / numel(z);
    maxdiff = max(abs(ZERO(:)));

    ret = perc == 0 && maxdiff < 10^(-perc);
    
%     for i = 1:N
%         ZERO = double(subs(z, s, prec*randn(size(s))));
%         if any(abs(ZERO(:)) > 10^(-prec))
%             ret = false;
%             return
%         end
%     end
%     ret = true;
end 

end