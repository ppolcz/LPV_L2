%% qLPV_ipend_MODEL_sim
%  
%  File: qLPV_ipend_MODEL_sim.m
%  Directory: workspace/3_ipend_qLPV
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2020. March 18. (2019b)
%

%%

m = 1;       % Mass of the rod [kg] 
M = 0.5;     % Mass of the car [kg]
L = 1;       % Length of the rod [m]
g = 10;      % Gravitational acceleration [m/s^2]
k = 0.5;     % Friction coefficient [??]

% Moment of inertia of a rod of length L and mass m, rotating about one
% end. This expression assumes that the rod is an infinitely thin (but
% rigid) wire.
I = m*L^2/3;

syms t y theta v omega real
x = [
    y
    v
    theta
    omega
    ];

p1 = sin(theta)/theta;
p2 = cos(theta);
p3 = omega*sin(theta);
sigma1 = (I+m*L^2)*(m+M) - m^2*L^2*p2^2;

A_qLPV = sigma1 \ [
    0 ,  sigma1      ,  0               ,  0 
    0 , -k*(I+m*L^2) , -m^2*L^2*g*p1*p2 , -m*L*(I+m*L^2)*p3
    0 ,  0           ,  0               ,  sigma1
    0 , -k*m*L*p2    , -(m+M)*m*g*L*p1  , -m^2*L^2*p2*p3
    ];

B_qLPV = sigma1 \ [
    0
    I+m*L^2
    0
    m*L*p2
    ];

f_sym = simplify(A_qLPV * x);
g_sym = B_qLPV;


f = matlabFunction(f_sym,'vars',{t x});
g = matlabFunction(g_sym,'vars',{t x});

%%

u = @(t) sin(4*t);
x0 = [0 0 0 0]';
[tt,xx] = ode45(@(t,x) f(t,x) + g(t,x)*u(t), [0 50], x0);

uu = u(tt);

ipend_simulate_Pi(tt,xx,uu)
