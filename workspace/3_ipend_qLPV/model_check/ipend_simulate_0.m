function ipend_simulate_0(t,x,u,r)
%% 
%
%  File: ipend_simulate_0.m
%  Directory: 4_gyujtemegy/11_CCS/Modellek/inverse_pendulum
%  Author: Peter Polcz (ppolcz@gmail.com) 
% 
%  Created on 2017.03.27. Monday, 13:46:12
%  Modified on 2019. November 14. (2019a)
%

dim = min(size(x));

% Frames per second
fps = 29;
Ts = 1/fps;

% Resampling x(t) using the interp1 function (`_rs` stands for `resampled`)
t_rs = t(1):Ts:t(end);
x_rs = interp1(t,x,t_rs);

% Handle input function
if nargin >= 3 
    if isnumeric(u) && numel(u) == numel(t)
        u_rs = interp1(t,u,t_rs);
    elseif isa(u,'function_handle')
        u_rs = u(t_rs);
        u = u(t);
    end
    input_title = 'Input';
else
    input_title = 'Input not given';
    u = t*0;
    u_rs = t_rs*0;
end

% Split the state vector
r_rs = x_rs(:,1);
if dim == 4
    phi_rs = x_rs(:,3);
    plot_1_data = x(:,1:2);
    plot_2_data = x(:,3:4);
elseif dim == 2
    phi_rs = x_rs(:,2);
    plot_1_data = x(:,1);
    plot_2_data = x(:,2);
else
    error 'The dimension of x must be 2 or 4'
end

% Handle reference input function
if nargin >= 4
    if isa(r,'function_handle')
        r = r(t);
    end
    plot_1_data = [plot_1_data r];
end

% Figure dimensions
l = 1;
cl = 1;    % length of chart
cw = 0.2;  % width of chart
rw = 0.05; % half width of the rod

% Sine and cosine functions of phi

s_rs = sin(phi_rs);
c_rs = cos(phi_rs);

% Coordinates of the free end of the rod
x1_rs = r_rs + l*s_rs;
y1_rs = l*c_rs;

% Corner points of the rod body
x11 = x1_rs + rw*c_rs;
y11 = y1_rs - rw*s_rs;

x12 = x1_rs - rw*c_rs;
y12 = y1_rs + rw*s_rs;

x21 = r_rs + rw*c_rs;
y21 = -rw*s_rs;

x22 = r_rs - rw*c_rs;
y22 = rw*s_rs;

% Plot helper for the chart
chartx = @(kk) [r_rs(kk)-cl/2, r_rs(kk)-cl/2, r_rs(kk)+cl/2, r_rs(kk)+cl/2, r_rs(kk)-cl/2];
charty = [-cw/2, cw/2, cw/2, -cw/2, -cw/2];

% Plot helper for the rod
rodx = @(kk) [x11(kk) x12(kk) x22(kk) x21(kk) x11(kk)];
rody = @(kk) [y11(kk) y12(kk) y22(kk) y21(kk) y11(kk)];

% Time display helper
simtime = @(kk) sprintf('t = %g [s]', round(t_rs(kk)));

% Dimension of the figure
fig = figure(19);

% [upper plots]
subplot(231), hold off;
plot(t, plot_1_data), xlim(t([1,end]))
title('Chart''s position and velocity'), grid on;

subplot(232), hold off;
plot(t, plot_2_data), xlim(t([1,end]))
title('Angle and angular velocity'), grid on;

subplot(233), hold off;
plot(t,u), xlim(t([1,end]))
title(input_title), grid on;

% [lower plot -- simulation movie]

% Delete previous plot (if exists)
delete(subplot(212));

% Create axis
ax = subplot(212); hold on

% Compute the optimal aspact ratio
figh = fig.Position(4) * ax.Position(4);
figw = fig.Position(3) * ax.Position(3);
plot_wph = figw/figh;
y_lim = [-2,2];
x_lim = max(y_lim * [-1;1] * plot_wph/2, 6) * [-1 1];

% Set coordinates limits
axis equal
xlim(x_lim);
ylim(y_lim);
grid on;

% Initialize plot
P_rod = fill(rodx(1),rody(1),[0.6350    0.0780    0.1840],'EdgeColor',[0.6350    0.0780    0.1840]);
P_chart = fill(chartx(1), charty,[0    0.4470    0.7410],'EdgeColor',[0    0.4470    0.7410]);
P_text = text(0,-1.5,simtime(1));
P_u = quiver(r_rs(1), -0.2, u_rs(1), 0, 'Color', [0.4660    0.6740    0.1880], 'LineWidth', 2);
drawnow

% Start simulation
tic;
for kk = 1:length(t_rs)
    P_rod.XData = rodx(kk);
    P_rod.YData = rody(kk);
    P_chart.XData = chartx(kk);
    P_u.XData = r_rs(kk);
    P_u.UData = u_rs(kk);
    P_text.String = simtime(kk);
    drawnow

    pause(Ts - toc)
    tic
end

end
