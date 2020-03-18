function [X,Y] = adapt_channels(X,Y)
%% adapt_channels
%  
%  File: adapt_channels.m
%  Directory: 1_PhD_projects/00_my_toolboxes/FinslerTools/v11/@PAffineMatrix
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2019. February 14.
%  Modified on 2019. November 05. (2019a)  [ EKKOR NEM MUKODOTT ]
%

%%

XisAffinMat = isa(X,'PAffineMatrix');
YisAffinMat = isa(Y,'PAffineMatrix');

if ~XisAffinMat && YisAffinMat
    [X,Y] = deal(Y,X);
    [XisAffinMat,YisAffinMat] = deal(YisAffinMat,XisAffinMat);
end

if XisAffinMat && ~YisAffinMat
    Y = PAffineMatrix(Y,1);
end

if numel(X.channels) ~= numel(Y.channels) || ~isempty(symvar(X.channels - Y.channels))

    U_channels = unique([ X.channels ; Y.channels ]);
    U_subsvars = unique([ X.subsvars ; Y.subsvars ]);

    if numel(U_channels) == numel(X.channels)
        U_channels = X.channels;
        U_subsvars = X.subsvars;
    end

    X = PAffineMatrix(X.get_channels__(U_channels),U_channels,U_subsvars);
    Y = PAffineMatrix(Y.get_channels__(U_channels),U_channels,U_subsvars);
end

end