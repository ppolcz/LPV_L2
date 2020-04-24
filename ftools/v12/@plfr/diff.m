function dLFR = diff(pLFR,p_cell,dp_cell,subsvars)
%%

dLFR = minlfr(pLFR*0);

for i = 1:numel(p_cell)
    try
        dLFR = dLFR + diff(pLFR.lfrtbx_obj,p_cell{i}.blk.names{1}) * dp_cell{i};
    catch e
        if strfind(e.message, 'does not exist in lfr-object')
            continue
        else
            getReport(e)
            break;
        end
    end
end

dp_plfr = plfr(vertcat(dp_cell{:}));

% if nargin < 4
%     subsvars = plfr.var_helper__(pLFR,dp_plfr);
% end
% 
% dLFR = dLFR.set_vars(subsvars(:));

end