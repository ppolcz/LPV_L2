function [ret] = pcz_get_plot_colors(N,I)
%% pcz_get_plot_colors
%  
%  File: pcz_get_plot_colors.m
%  Directory: 7_ftools/utilities/plotting_tools
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2019. June 20. (2019a)
%

if nargin <= 0
    N = 10;
end

DO_PLOT = true;
if isempty(N)
    N = 10;
    DO_PLOT = false;
end

if nargin <= 1
    I = 1:N;
end

%%

fig = figure;
P = plot(0:1,repmat(1:N,[2,1]), 'LineWidth',5);

ret = vertcat( P.Color );

ret = ret(I,:);

if ~DO_PLOT
    delete(fig)
end

end