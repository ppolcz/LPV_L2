function constKerT = P_constKerT_for_LFR(lfr)
%% P_constKerT_for_LFR
%  
%  File: P_constKerT_for_LFR.m
%  Directory: 1_PhD_projects/00_my_toolboxes/FinslerTools/v11
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2019. March 10.
%

%%


PI_lfr = lfr;
if isa(PI_lfr,'lfr')
    PI_lfr = plfr(lfr);
end


N = PI_lfr.np^2 + PI_lfr.ny - PI_lfr.nu;

samples = randn(PI_lfr.np,N);

samples_cell = num2cell(samples,1);

behely_cell = cellfun(@(p_num) { PI_lfr(p_num)' }, samples_cell);

M = vertcat(behely_cell{:});

constKerT = null(M)';


end