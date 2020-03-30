function [ret] = helper_check_stability(A_fh, B_fh, C_fh, D_fh, p_lim)
% 
%  File: helper_check_stability.m
%  Directory: 8_published/LPV_L2/workspace
%  Author: Peter Polcz (ppolcz@gmail.com) 
% 
%  Created on 2018. November 19.
%  Major review on 2020. March 26. (2019b)
%

%%

if nargin == 5
    [A_fh, B_fh, C_fh, D_fh] = helper_fh2fh(A_fh, B_fh, C_fh, D_fh);
else
    p_lim = B_fh;
end

pcz_dispFunctionSeparator
pcz_dispFunction('Stability check in random points')
ok = 1;
for i = 1:100
    p_num = rand(size(p_lim,1),1).*(p_lim(:,2)-p_lim(:,1))+p_lim(:,1);
    p_cell = num2cell(p_num);
    if pcz_failif(all(real(eig(A_fh(p_cell{:}))) < 0), 'A(p) is NOT stable in %s', pcz_num2str(p_num))
        ok = 0;
        eig(A_fh(p_cell{:}))
    end
    
    % sys = ss(A_fh(p_cell{:}), B_fh(p_cell{:}), C_fh(p_cell{:}), D_fh(p_cell{:}));
    % pcz_dispFunction2('If p = [%6.3g,%6.3g,%6.3g], norm = %g.', p_num, norm(sys,Inf))
end

P_v = P_ndnorms_of_X(p_lim);
pcz_dispFunction('Stability check in corner points')
for i = 1:size(P_v,1)
    p_num = P_v(i,:)';
    if pcz_failif(all(real(eig(A_fh(p_cell{:}))) < 0), 'A(p) is NOT stable in %s', pcz_num2str(p_num))
        ok = 0;
    end
end

if ok
    pcz_info(1,'A(p) was stable in all sample points')
end

end