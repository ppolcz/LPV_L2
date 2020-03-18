function ipend_simulate_Pi(t,x,varargin)
%% 
%  
%  File: ipend_simulate_Pi.m
%  Directory: 4_gyujtemegy/11_CCS/Modellek/inverse_pendulum
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2017.03.27. Monday, 14:53:36
%  Modified on 2019. November 15. (2019a)
%

dim = min(size(x));

if dim == 4
    x(:,3) = x(:,3) + pi;
elseif dim == 2
    x(:,2) = x(:,2) + pi;
end    
ipend_simulate_0(t,x,varargin{:})

end