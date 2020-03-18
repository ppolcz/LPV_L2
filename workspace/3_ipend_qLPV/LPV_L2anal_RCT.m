function LPV_L2anal_RCT(A_fh, B_fh, C, D, K, p_lim)
%% LPV_L2anal_RCT
%
%  File: LPV_L2anal_RCT.m
%  Directory: 1_PhD_projects/22_Hinf_norm/LPV_TAC_v3_newmod
%  Author: Peter Polcz (ppolcz@gmail.com)
%
%  Created on 2018. November 19.
%


%%
% Automatically generated stuff

global SCOPE_DEPTH VERBOSE LATEX_EQNR
SCOPE_DEPTH = 0;
VERBOSE = 1;
LATEX_EQNR = 0;

TMP_lfqlHdooTxNFOfYScyqG = pcz_dispFunctionName('Robust Control Toolbox');

[nz,nx] = size(C);
[~,nw] = size(D);
np = size(p_lim,1);

p_cell = cell(1,np);
for i = 1:np
    name = [ 'p' num2str(i) ];
    p_nom = sum(p_lim(i,:))/2;
    p_cell{i} = ureal(name, p_nom, 'Range', p_lim(i,:));
end

A_unc = A_fh(p_cell{:});
B_unc = B_fh(p_cell{:});
C_unc = umat(C);
D_unc = umat(D);

Ak_unc = A_unc - B_unc*K;

sys_unc = ss(Ak_unc,B_unc,C_unc,D_unc);
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

    bounds = [wcg.LowerBound , wcg.UpperBound];

    pcz_dispFunction('Solver time: <strong>%g</strong>', toc(TMP_KZWeXiYFmdpQdsgidKeG))
    pcz_dispFunction(2, 'Bounds: [<strong>%g</strong>,%g] = [%g,%g] dB ', bounds, 20*log10(bounds))
    parnames = fields(wcu);
    for fld = 1:numel(parnames)
        pcz_dispFunction('Worst case uncertainty: %s = %g', parnames{fld}, wcu.(parnames{fld}))
    end
    pcz_dispFunction('Critical frequency: %g', wcg.CriticalFrequency)
    pcz_dispFunctionEnd(TMP_KZWeXiYFmdpQdsgidKeG);

    % -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -

%     figure('Position', [  520.000 ,  492.000 , 1119.000 ,  359.000 ], 'Color', [1 1 1])
%     opts = wcOptions('VaryFrequency','on','Display','on');
%     [wcg,~,info] = wcgain(sys_unc,linspace(0.1,10,10),opts)
%     subplot(121), semilogx(info.Frequency,info.Bounds), grid on
% 
%     title('Worst-Case Gain vs. Frequency')
%     ylabel('Gain')
%     xlabel('Frequency')
%     legend('Lower bound','Upper bound','Location','northwest')
% 
%     wc_sys = usubs(sys_unc,wcu);
%     subplot(122), sigma(wc_sys), grid on

% -------------------------------------------------------------------------

%{

wcmargin(sys_unc)

%}


%{
% -------------------------------------------------------------------------
% Compute analytically the DC-gain

    [A,C,B,D] = pcz_split_matrix(F_fh(wcu.p1),[nx,nz],[nx,nw], 'RowWise', false);
    wc_sys_ = tf(ss(A,B,C,D));

    pcz_dispFunction('DC gain: %g', dcgain(wc_sys_))

% -------------------------------------------------------------------------
%}


pcz_dispFunctionEnd(TMP_lfqlHdooTxNFOfYScyqG)

end


%%

% sys = ss(A_fh(wcu.p1,wcu.p2,wcu.p3)-B_fh(wcu.p1,wcu.p2,wcu.p3)*K, B_fh(wcu.p1,wcu.p2,wcu.p3), C, D);
% bode(sys)
% grid on
