function [X,Y] = adapt_channels(X,Y)
%% adapt_channels
%  
%  File: adapt_channels.m
%  Directory: 1_PhD_projects/00_my_toolboxes/FinslerTools/v11/@PGenAffineMatrix
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2019. February 14.
%

%%

XisAffinMat = isa(X,'PGenAffineMatrix');
YisAffinMat = isa(Y,'PGenAffineMatrix');

if ~XisAffinMat && YisAffinMat
    [X,Y] = deal(Y,X);
    [XisAffinMat,YisAffinMat] = deal(YisAffinMat,XisAffinMat);
end

if XisAffinMat && ~YisAffinMat
    Y = PGenAffineMatrix(Y,1);
end

if numel(X.channels) ~= numel(Y.channels) || ~isempty(symvar(X.channels - Y.channels))

    U_channels = unique([ X.channels ; Y.channels ]);
    U_subsvars = unique([ X.subsvars ; Y.subsvars ]);

    if numel(U_channels) == numel(X.channels)
        U_channels = X.channels;
        U_subsvars = X.subsvars;
    end

    X = PGenAffineMatrix(X.get_channels__(U_channels),U_channels,U_subsvars);
    Y = PGenAffineMatrix(Y.get_channels__(U_channels),U_channels,U_subsvars);
end

end