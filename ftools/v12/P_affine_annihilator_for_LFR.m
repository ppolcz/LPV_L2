function [N_lfr,samples,sampleserr,err] = P_affine_annihilator_for_LFR(PI, vars_lfr, varargin)
%% P_affine_annihilator_for_LFR
%  
%  File: P_affine_annihilator_for_LFR.m
%  Directory: 1_PhD_projects/00_my_toolboxes/FinslerTools/v11
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2019. March 10.
%

TMP_fzhqgrIAPDdmVUZqZFHF = pcz_dispFunctionName;

%%

if isa(PI,'lfr')
    PI = plfr(PI);
end

args.tol = 1e-7;
args.proj = [];
args.symbolic = 0;
args.precision = 10;
args.samples = [];
args.sampleserr = [];
args.lims_or_np = PI.np;
args = parsepropval(args,varargin{:});

%%

% 2019.11.21. (november 21, csütörtök), 12:04
vars = sym(plfr(vars_lfr));

m = PI.ny;
np = size(vars_lfr,1);

ABn = eye(m*(np+1));

An = ABn(:,1:m);
Bn = ABn(:,m+1:end);
Cn = kron(ones(np,1),eye(m));
Dn = zeros(m*np);

blk = vars_lfr.blk;
blk.desc(1:2,:) = ones(2,np)*m;

N_kiemelt_lfr = lfr(Dn,Cn,Bn,An,blk);

N_times_PI = N_kiemelt_lfr*PI.lfrtbx_obj;

% MODOSITVA: 2019.11.21. (november 21, csütörtök), 12:06
% `N_times_PI' --> plfr(N_times_PI',vars_lfr)'
[N_ThetaT,~,samples,sampleserr,err] = P_constKer_for_GSS(plfr(N_times_PI',vars), 'proj', args.proj,...
    'samples', args.samples, 'sampleserr', args.sampleserr, 'lims', args.lims_or_np, 'tol', args.tol);

N_Theta = N_ThetaT';

if args.symbolic == 1
    % beautifying
    N_Theta = rref(N_Theta);

    % !! rounding !!
    N_Theta = double(pcz_find_recdec(N_Theta, 'tolerance', 10^(-args.precision)));
end

N_lfr = plfr(N_Theta * N_kiemelt_lfr, vars);

pcz_dispFunctionEnd(TMP_fzhqgrIAPDdmVUZqZFHF);

end



% blk = struct;
% blk.desc = [
%     ones(2,np)*m  % dimension
%     ones(5,np)    % minden csupa 1
%     ones(1,np)*2  % gondolom a bound tipus (minmax)
%     -ones(1,np)   % lower bound
%     ones(1,np)    % upper bound
%     zeros(1,np)   % nominal value
%     ];
% blk.names = cellfun(@(p) {char(p)}, num2cell(vars_lfr)).';
% 
% indices = reshape(1:(np*m+m),[m,np+1])';
% N_num = PAffineMatrix(N_Theta(:,indices(:)),[1;vars_lfr]);
% 
% N_sym = P_affine_annihilator(Pi,xp,vars_lfr,'sym',1);
% N_sym = N_sym.set_subsvars(vars_lfr);
% 
% 
% PI_fh = matlabFunction(PI, 'vars', {vars_lfr});
% pcz_fhzero_report(@(p) N_sym(p)*PI_fh(p),vars_lfr,1e-10,100)
% pcz_fhzero_report(@(p) N_num(p)*PI_fh(p),vars_lfr,1e-10,100)
% 
