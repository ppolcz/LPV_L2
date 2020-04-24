%% LPV4D_main
%
%  File: LPV4D_main.m
%  Directory: workspace/1_comp_LPV
%  Author: Peter Polcz (ppolcz@gmail.com)
%
%  Created on 2020. March 20. (2019b)
%

G_reset
P_init(12)

setenv('RUN_ID', num2str(pcz_runID(mfilename)))
logger = Logger(['results/' mfilename '-output.txt']);
TMP_vcUXzzrUtfOumvfgWXDd = pcz_dispFunctionName;
pcz_dispFunction2('Run ID = %s', getenv('RUN_ID'));


%%

P_generate_symvars_v10(4,3,2,2)

A_fh = @(p1,p2,p3) [
    -3+p1       3+p1/(p3^2 + 0.5*p1 + 1)     0.1*p3   p3^2*p2
    0          -1-p2^2                       5        0
    -1/(5-p2)  0                             -4+p1    0
    0          0.1                           0        -5+1/(p1 + 2)
    ];

C_fh = @(p1,p2,p3) [
    1/(5-p2)   0                             0        0
    0          0                             p1+1     0
    ];

B_fh = @(p1,p2,p3) [
    0      0
    1+p2^2 0
    0      0
    0      2+p1/(p1 + 2)
    ];

D_fh = [
    0      0
    0      0
    ];

p_lim = [
    -1 2
    -1 2
    0 2
    ];

% LPV_quick_check_stability(A_fh, B_fh, C_fh, D, p_lim)

% Basis functions for the grid-based methods
bases = [
                  1
      monomials(p,1) % Requires SOS Tools
      monomials(p,2) % Requires SOS Tools
      monomials(p,3) % Requires SOS Tools
            1/(5-p2)
           p1/(p1+2)
  p1/(p3^2+0.5*p1+1)
    ];
bases_Jac = jacobian(bases,p);


%% N logarithmically equidistant points in a decade
for Scale_dp_lim = setdiff(unique([
%         logspace(0,1,11)'
%         logspace(1,2,11)'
%         ...
        logspace(0,1,5)'
        logspace(1,2,5)'
        logspace(2,3,5)'
%         ...
%         logspace(3,4,4)'
%         logspace(4,5,4)'
%         ...
%         logspace(4,5,2)'
        ])',[])
%%

modelname = sprintf('model1_LPV4D_x%g', Scale_dp_lim);

dp_lim = [
    -10 10
    -1 1
    -5 5
    ]*Scale_dp_lim;


p_lims_comp = p_lim;

pdp_lims_comp = [
    p_lims_comp
    dp_lim
    ];

%%

% method1_RCT(modelname,A_fh,B_fh,C_fh,D_fh,p_lim)

% method0_grid_LPVTools(modelname,A_fh,B_fh,C_fh,D_fh,p_lim,dp_lim,bases,bases_Jac,[5 5 5]);
%
% % Greedy grid
% method0_grid_Wu1995(modelname,A_fh,B_fh,C_fh,D_fh,p_lim,dp_lim,bases,bases_Jac,'res_max',5);
method0_grid_Wu1995(modelname,A_fh,B_fh,C_fh,D_fh,p_lim,dp_lim,bases,bases_Jac,'res_max',13);
% method0_grid_Wu1995(modelname,A_fh,B_fh,C_fh,D_fh,p_lim,dp_lim,bases,bases_Jac,'res_max',15);
%
% % As proposed by Wu (1995,1996)
% method0_grid_Wu1995(modelname,A_fh,B_fh,C_fh,D_fh,p_lim,dp_lim,bases,bases_Jac,'res_max',5,'delta',1e-6,'T',1e6);
% method0_grid_Wu1995(modelname,A_fh,B_fh,C_fh,D_fh,p_lim,dp_lim,bases,bases_Jac,'res_max',5,'delta',1e-5,'T',100000);
% method0_grid_Wu1995(modelname,A_fh,B_fh,C_fh,D_fh,p_lim,dp_lim,bases,bases_Jac,'res_max',5,'delta',1e-4,'T',10000);
% method0_grid_Wu1995(modelname,A_fh,B_fh,C_fh,D_fh,p_lim,dp_lim,bases,bases_Jac,'res_max',5,'delta',1e-3,'T',1000);
% method0_grid_Wu1995(modelname,A_fh,B_fh,C_fh,D_fh,p_lim,dp_lim,bases,bases_Jac,'res_max',5,'delta',1e-2,'T',100);
% method0_grid_Wu1995(modelname,A_fh,B_fh,C_fh,D_fh,p_lim,dp_lim,bases,bases_Jac,'res_max',5,'delta',1e-1,'T',100);
% method0_grid_Wu1995(modelname,A_fh,B_fh,C_fh,D_fh,p_lim,dp_lim,bases,bases_Jac,'res_max',5,'delta',1e-2,'T',10);
% method0_grid_Wu1995(modelname,A_fh,B_fh,C_fh,D_fh,p_lim,dp_lim,bases,bases_Jac,'res_max',5,'delta',1e-1,'T',10);
%
% method2_descriptor_primal(modelname,A_fh,B_fh,C_fh,D_fh,p_lim,dp_lim)
% method2_descriptor_dual(modelname,A_fh,B_fh,C_fh,D_fh,p_lim,dp_lim,1)
% method2_descriptor_dual(modelname,A_fh,B_fh,C_fh,D_fh,p_lim,dp_lim,0)
%
% % IQC/LFT approaches for LPV with rate-bounded parameters
% method3_IQC_LFT_LPVMAD(modelname,A_fh,B_fh,C_fh,D_fh,p_lim,dp_lim);
% method3_IQC_LFT_LPVTools(modelname,A_fh,B_fh,C_fh,D_fh,p_lim,dp_lim);
%
% method4_authors_old_symbolical(modelname,A_fh,B_fh,C_fh,D_fh,p_lim,dp_lim,...
%     'minlfr', true);
%
% % Imported variables to the base workspace: Q, dQ, PI_x, gamma
% method5_proposed_approach(modelname,A_fh,B_fh,C_fh,D_fh,p_lim,dp_lim,p_lims_comp,pdp_lims_comp,...
%     'minlfr', true);

end

pcz_dispFunctionEnd(TMP_vcUXzzrUtfOumvfgWXDd);
logger.stoplog

return


%% Plot overall results

Entries = {
    '$\hat\gamma$, grid $5\times5\times5$ $V(x,p)$' '^--'  5 true
    '$\hat\gamma$, grid $7\times7\times7$ $V(x,p)$' '^--'  5 false
    '$\hat\gamma$, grid $13\times13\times13$ $V(x,p)$' '^--'  5 true
    'upper bound $\gamma$, descriptor'       '.--'  30 true
    'upper bound $\gamma$, LPVMAD'           '.--'  40 true
    'upper bound $\gamma$, lpvwcgain'        '.--'  20 true
    'upper bound $\gamma$, Finsler (old)'    's--'  10 true
    'upper bound $\gamma$, Finsler (new)'    '.k--' 10 true
    };


Res = [
1	1.26	1.59	1.78	2	2.15	2.51	3.16	3.98	4.64	5.01	5.6	6.31	7.94	10	17.8	21.5	31.6	46.4	56.2	100	177.8	215.4	316.2	464.2	562.3	1000	1.00E+04	1.00E+05
2.03523957595921	NaN	NaN	2.03523957668091	NaN	NaN	NaN	2.05545642480682	NaN	NaN	NaN	2.14528742985676	NaN	NaN	2.26678243850142	2.3896684602428	NaN	2.49177817818668	NaN	2.56546022948363	NaN	2.64263560252844	NaN	2.65992389825714	NaN	2.6699170221645	2.67562535952935	NaN	NaN
2.03523957600935	NaN	NaN	2.0375360506381	NaN	NaN	NaN	2.07836868819199	NaN	NaN	NaN	2.170878918141	NaN	NaN	2.30270111649155	2.44264158416503	NaN	2.56064016506232	NaN	2.64514028230515	2.69983381045177	2.7331560617445	NaN	2.75275034058177	NaN	2.76405264863283	2.77049996025163	NaN	NaN
2.0352395759488	NaN	NaN	2.03877065378142	NaN	NaN	NaN	2.07952106297757	NaN	NaN	NaN	2.17289013401731	NaN	NaN	2.30591180984638	2.4454328668872	NaN	2.56234109618297	NaN	2.64603162112553	2.70038247672717	2.73350365183496	NaN	2.75300025125219	NaN	2.76421448395179	2.7705980545329	NaN	NaN
2.10564666849556	NaN	NaN	2.20063554509286	NaN	NaN	NaN	2.36050330073879	NaN	NaN	NaN	2.55228514512871	NaN	NaN	2.69285796471826	2.75261587698898	NaN	2.7723471343514	NaN	2.77834860730295	2.78009235130869	2.78058262438168	NaN	2.7807241834886	NaN	2.78075685811572	2.7807585761487	2.78075936139437	2.78075911067564
2.24860849805858	2.33489739789273	2.46624325313265	2.54833429486015	2.63818392553585	NaN	2.80114519282558	2.81918562097354	2.81977241055676	NaN	2.81860414924753	2.81794253677084	2.81764609080954	2.82018392520324	2.82033877316899	2.820254204693	NaN	2.81888218089099	NaN	2.81994196725838	2.82026959237618	2.82072910307767	NaN	2.82318853468453	NaN	2.82031935321688	2.81963134967875	2.82739191302986	2.87240145997068
2.36441711121937	2.4367725457923	2.53209068828656	2.59231592853063	2.65813032943977	NaN	2.79533678925321	2.81854965249846	2.81852067142831	NaN	2.81839227737	2.81812896808011	2.819047448031	2.81921458920945	2.81884398630243	2.81819053246449	NaN	2.81913762542185	NaN	2.81917041397138	2.81924566906359	2.81932602818647	NaN	2.81930701683456	NaN	2.81813077694576	2.81813194525962	2.81862859240125	2.81851428506674
2.03524484995885	NaN	NaN	2.05119698597003	NaN	NaN	NaN	2.11778843658819	NaN	NaN	NaN	2.23185778909705	NaN	NaN	2.37246558356987	2.50699805713682	NaN	2.61180781721747	NaN	2.68530090387105	2.72932634897617	2.75467733095872	NaN	2.76579438514387	NaN	2.77425576002343	2.77458783128181	NaN	NaN
2.03524102481956	NaN	NaN	2.05120659809694	NaN	NaN	NaN	2.11755668912915	NaN	NaN	NaN	2.23104291531253	NaN	NaN	2.37250963593153	2.50689326637037	NaN	2.61089501020792	NaN	2.68469457975552	2.7291213977834	2.75319928434571	NaN	2.76510884284798	NaN	2.77316765539124	2.7768342796527	2.694649013175	NaN
    ];

s = Res(1,:)';
data = num2cell(Res(2:end,:)',1);

HARD_LOWER = norm(ss(A_fh(2,2,2),B_fh(2,2,2),C_fh(2,2,2),D_fh),Inf);

figure(1), delete(gca), hold on
for i = 1:numel(data)
    I = ~isnan(data{i});
    if Entries{i,4}
        plot(s(I), data{i}(I),Entries{i,2},'MarkerSize',Entries{i,3})
    end
end
plot(s,s*0+HARD_LOWER,'r','LineWidth',4);
set(gca,'xscale','log')
Leg = legend(Entries{[Entries{:,4}],1});
Leg.Location = 'southeast';
Leg.Interpreter = 'latex';
Leg.FontSize = 14;

grid on

xlim([1,1000])
ylim([2,2.9])

Logger.latexify_axis(14)
Logger.latexified_labels(gca,16,'Value of $\alpha$','Computed $\gamma$')

%{

print('results_stored/model1_LPV4D-plot-3.pdf','-dpdf')

%}