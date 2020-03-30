function method1_RCT(modelname, A_fh, B_fh, C_fh, D_fh, p_lim)
%%
%  File: method1_grid_RCT.m
%  Directory: 8_published/LPV_L2/workspace
%  Author: Peter Polcz (ppolcz@gmail.com) 
% 
%  Created on 2018. November 19.
%  Major review on 2020. March 26. (2019b)


%%
% Automatically generated stuff

TMP_lfqlHdooTxNFOfYScyqG = pcz_dispFunctionName('Robust Control Toolbox');

np = size(p_lim,1);
p_nom = sum(p_lim,2)/2;

p_cell = cell(1,np);
for i = 1:np
    name = [ 'p' num2str(i) ];
    p_cell{i} = ureal(name, p_nom(i), 'Range', p_lim(i,:));
end

[A_unc,B_unc,C_unc,D_unc] = helper_convert(A_fh,B_fh,C_fh,D_fh,p_cell,'umat');

sys_unc = ss(A_unc,B_unc,C_unc,D_unc);
simplify(sys_unc,'full');

% -------------------------------------------------------------------------
% Bode plot

    figure
    bopts = bodeoptions;
    bopts.MagUnits = 'abs';
    bodeplot(gridureal(sys_unc,30), bopts), grid on

% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% Peak gain

    Nr_samples = 50;
    [peak_gain,freki] = getPeakGain(gridureal(sys_unc,Nr_samples));

    [peak_gain,I] = max(peak_gain);
    freki(I);

% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% Worst case gain

    TMP_KZWeXiYFmdpQdsgidKeG = pcz_dispFunctionName('Pure wcgain');
    [wcg,wcu,info] = wcgain(sys_unc);

    pcz_dispFunction2(evalc('display(info)'));
    
    bounds = [wcg.LowerBound , wcg.UpperBound];

    Overall_Time = toc(TMP_KZWeXiYFmdpQdsgidKeG);
    
    pcz_dispFunction('Solver time: <strong>%g</strong>', Overall_Time)
    pcz_dispFunction(2, 'Bounds: [<strong>%g</strong>,%g] = [%g,%g] dB ', bounds, 20*log10(bounds))
    parnames = fields(wcu);
    for fld = 1:numel(parnames)
        pcz_dispFunction('Worst case uncertainty: %s = %g', parnames{fld}, wcu.(parnames{fld}))
    end
    pcz_dispFunction('Critical frequency: %g', wcg.CriticalFrequency)
    pcz_dispFunctionEnd(TMP_KZWeXiYFmdpQdsgidKeG);

    
    store_results('Results_All.csv', modelname, bounds(1), bounds(2), 0, Overall_Time, ...
        '[no info]', 'Robust Control Toolbox')
    
% -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -

pcz_dispFunctionEnd(TMP_lfqlHdooTxNFOfYScyqG)

end
