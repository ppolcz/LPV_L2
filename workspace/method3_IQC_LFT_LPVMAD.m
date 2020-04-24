function method3_IQC_LFT_LPVMAD(modelname, A_fh, B_fh, C_fh, D_fh, p_lim, dp_lim)
%% LPV_L2anal_IQCToolbox
%  
%  File: LPV_L2anal_IQCToolbox.m
%  Directory: 1_PhD_projects/22_Hinf_norm/qLPV_ipend
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2020. March 13. (2019b)
%

global LPVMAD_vars
LPVMAD_vars.Solver_Time = 0;

%%

TMP_quNJgGJNllaEMSewPAvy = pcz_dispFunctionName('IQCToolbox');

% Translate the parameter bounds such that the origin be in the middle.
p_nom = sum(p_lim,2)/2;
p_centered_lim = p_lim - p_nom;

% Generate translated lfr variables
[~,p_cell] = pcz_generateLFRStateVector('p',p_centered_lim);
p_centered_cell = cellfun(@(p,p0) { p + p0 }, p_cell(:), num2cell(p_nom));

[~,~,~,~,~,~,M_lfr,nx,~,nu,ny] = helper_convert(A_fh,B_fh,C_fh,D_fh,p_centered_cell);

M_lfr = minlfr(M_lfr);

Fij = [
    M_lfr.d M_lfr.c
    M_lfr.b M_lfr.a
    ];

F = cell(3);
[F{:}] = pcz_split_matrix(Fij,[nx,ny,Inf],[nx,nu,Inf],'RowWise',0);

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

A = F{1,1} - eye(size(F{1,1}))*eps;

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

[iqc,~] = iqcuc(iqc,'ltvs',M_lfr.blk.desc(1,:),[p_centered_lim' dp_lim'],'box');
 
nr_blocks = numel(M_lfr.blk.desc(1,:));
Pole = -1 * ones(1,nr_blocks);
Length = 2 * ones(1,nr_blocks);

iqc = set(iqc,1,'Pole',Pole);
iqc = set(iqc,1,'Length',Length);
iqc = set(iqc,1,'Relax',1);
iqc = set(iqc,1,'Active',1);

tic
iqc = iqcsolve(iqc);
Overall_Time = toc;
toc

get(iqc)

gamma = get(iqc,'gopt');
if isempty(gamma)
    gamma = 0;
end

pcz_dispFunction(2,'Model: %s', modelname);
pcz_dispFunction(2,'Overall time: %g', Overall_Time);
pcz_dispFunction(2,'<strong>Gamma = %g</strong> ', gamma);
pcz_dispFunction(' ');

info = [ 'Pole: [' num2str(Pole) '], Length: [' num2str(Length) '], Relax: ' num2str(get(iqc,1,'Relax')) ', Active: ' num2str(get(iqc,1,'Active')) ];

store_results('Results_All',modelname,0,gamma,LPVMAD_vars.Solver_Time,Overall_Time,info,'LPVMAD - ltvs')

%%
pcz_dispFunctionEnd(TMP_quNJgGJNllaEMSewPAvy);

end