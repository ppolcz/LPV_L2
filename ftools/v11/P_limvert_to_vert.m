function [P_v, P_Nr, Mode_str] = P_limvert_to_vert(plim, nrvars)
%% P_limvert_to_vert
%  
%  File: P_limvert_to_vert.m
%  Directory: 1_PhD_projects/00_my_toolboxes/FinslerTools/v11
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2019. February 15.
%

%%

Mode_plim = 0;

if size(plim,2) == 2 && all(plim * [-1 ; 1] > 0) && nrvars == size(plim,1)
    [P_v,~,~] = P_ndnorms_of_X(plim);
    Mode_plim = 1;
else
    P_v = plim;
end

P_Nr = size(P_v,1);


Mode_str = { 'vertices are given' 'rectangular region: limits are given' };
Mode_str = Mode_str{Mode_plim + 1};


end