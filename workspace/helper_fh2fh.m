function [A_fh,B_fh,C_fh,D_fh,AC_fh,BD_fh,M_fh] = helper_fh2fh(A_fh, B_fh, C_fh, D_fh)
%% helper_fh2fh
%
%  File: helper_fh2fh.m
%  Directory: 8_published/LPV_L2/workspace
%  Author: Peter Polcz (ppolcz@gmail.com)
%
%  Created on 2020. March 25. (2019b)
%

%%

if isnumeric(A_fh)
    A_fh = @(varargin) A_fh;
end

if isnumeric(B_fh)
    B_fh = @(varargin) B_fh;
end

if isnumeric(C_fh)
    C_fh = @(varargin) C_fh;
end

if isnumeric(D_fh)
    D_fh = @(varargin) D_fh;
end

AC_fh = @(varargin) [
    A_fh(varargin{:})
    C_fh(varargin{:})
    ];

BD_fh = @(varargin) [
    B_fh(varargin{:})
    D_fh(varargin{:})
    ];

M_fh = @(varargin) [ AC_fh(varargin{:}) BD_fh(varargin{:}) ];

end