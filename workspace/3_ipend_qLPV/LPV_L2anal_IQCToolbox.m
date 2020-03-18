function LPV_L2anal_IQCToolbox(A_fh, B_fh, C, D, K, x_lim, p_lim, dp_lim)
%% LPV_L2anal_IQCToolbox
%  
%  File: LPV_L2anal_IQCToolbox.m
%  Directory: 1_PhD_projects/22_Hinf_norm/qLPV_ipend
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2020. March 13. (2019b)
%

%%
% Automatically generated stuff

global SCOPE_DEPTH VERBOSE LATEX_EQNR 
SCOPE_DEPTH = 0;
VERBOSE = 1;
LATEX_EQNR = 0;

try c = evalin('caller','persist'); catch; c = []; end
persist = Persist(mfilename('fullpath'), c); clear c; 
persist.backup();
%clear persist

%%

% Find out np, nw and nz
[nz,nx] = size(C);
[~,nw] = size(D);
np = size(p_lim,1);

% Generate symbolic lfr variables
[p_lfr,p_cell] = pcz_generateLFRStateVector('p',p_lim);
[dp_lfr,dp_cell] = pcz_generateLFRStateVector('dp',dp_lim);


p_cell{1} = p_cell{1} + 1;
p_cell{2} = p_cell{2} + 1;
p_lim(1:2,1:2) = p_lim(1:2,1:2) - 1;


% LFR realization of matrix B(p):
B_lfr = B_fh(p_cell{:});

% LFR realization of matrix A(p) - B(p)*K:
Ak_lfr = A_fh(p_cell{:}) - B_lfr*K;

F_lfr_initial = [
    Ak_lfr B_lfr
    C      D
    ];

F_lfr = minlfr(F_lfr_initial);
% P_lfr = P_lfrdata_v6(F_lfr_initial)

%%

m = size(F_lfr.a,1);

Fij = [
    F_lfr.d F_lfr.c
    F_lfr.b F_lfr.a
    ];
F = cell(3);

F{1,1} = Fij(1:nx       ,1:nx); F{1,2} = Fij(1:nx       ,nx+1:nx+nw); F{1,3} = Fij(1:nx       ,nx+nw+1:end);
F{2,1} = Fij(nx+1:nx+nw ,1:nx); F{2,2} = Fij(nx+1:nx+nw ,nx+1:nx+nw); F{2,3} = Fij(nx+1:nx+nw ,nx+nw+1:end);
F{3,1} = Fij(nx+nw+1:end,1:nx); F{3,2} = Fij(nx+nw+1:end,nx+1:nx+nw); F{3,3} = Fij(nx+nw+1:end,nx+nw+1:end);

% A possible shorthand:
% [F{:}] = pcz_split_matrix(Fij, [nx nz m], [nx nw m], 'RowWise', false);


% In my model
%
% dx = F11 x + F12 w + F13 Pi
%  z = F21 x + F22 w + F23 Pi
%  y = F31 x + F32 w + F33 Pi

% In the IQC Toolbox model [THIS ONE USED IN THIS CASE]
%
% dx = F11 x + F13 Pi + F12 w
%  y = F31 x + F33 Pi + F32 w
%  z = F21 x + F23 Pi + F22 w

A = F{1,1} - eye(3)*eps;

B = [
    F{1,3} F{1,2}
    ];

C = [
    F{3,1}
    F{2,1}
    ];

D = [
    F{3,3} F{3,2}
    F{2,3} F{2,2}
    ];

M = ss(A,B,C,D);

iqc = iqcpb(M);

[iqc,ucpointer3] = iqcuc(iqc,'ltvs',F_lfr.blk.desc(1,:),[p_lim' dp_lim'],'box');

% Best results for:
% iqc = set(iqc,1,'Pole',[-2 -2 -2]);
% iqc = set(iqc,1,'Length',[2 2 2]);

 
iqc = set(iqc,1,'Pole',[-2 -2 -2]);
iqc = set(iqc,1,'Length',[2 2 2]);
% iqc = set(iqc,1,'Length',[1 1 1]);
iqc = set(iqc,1,'Relax',1);
iqc = set(iqc,1,'Active',1);

tic
iqc = iqcsolve(iqc);
IQC_solution_time = toc;
toc

get(iqc)

gamma = get(iqc,'gopt');
if isempty(gamma)
    gamma = 0;
end

info = [ 'Pole: [' num2str(get(iqc,1,'Pole')) '], Length: [' num2str(get(iqc,1,'Pole')) '], Relax: ' num2str(get(iqc,1,'Relax')) ', Active: ' num2str(get(iqc,1,'Active')) ];

%% Store results

s.Stamp = persist.stamp;
s.x1_max = x_lim(1,2);
s.x2_max = x_lim(2,2);
s.x3_max = x_lim(3,2);
s.Gamma = gamma;
s.Solver_Time = IQC_solution_time;
s.Solver_info = info;
s.Method = 'IQCToolbox - ltvs';
s.Use_MinLFR = 0;

Results_csv = persist.file_simple('txt','IQCToolbox_results.csv');
pcz_dispFunction('IQCToolbox_results.csv: `%s''', Results_csv);

if exist(Results_csv, 'file')
    Results = readtable(Results_csv,'Delimiter','|');
end

if ~exist('Results', 'var')
    Cols = { 'Stamp', 'x1_max', 'x2_max', 'x3_max', 'Gamma', 'Solver_Time', 'Solver_info', 'Method', 'Use_MinLFR' };
    Results = cell2table(cell(0,numel(Cols)), 'VariableNames', Cols);
end

s

Results
struct2table(s)

Results = [ Results ; struct2table(s) ];

% save(Results_mat,'Results');
writetable(Results,Results_csv,'Delimiter','|');
clear Results


%%
persist.stoplog;

end