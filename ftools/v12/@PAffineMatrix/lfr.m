function syslfr = lfr(N,pl_cell)
%% lfr
%  
%  File: lfr.m
%  Directory: 7_ftools/ftools/v11/@PAffineMatrix
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2019. November 05. (2019a)
%  Major review on 2020. March 30. (2019b)
%

%%

if nargin > 1 && iscell(pl_cell) && isa(pl_cell{1},'lfr')
    
    % 2020.03.30. (március 30, hétfő), 22:56
    paffmat_names = cellfun(@(s) { char(s) }, num2cell(N.channels(:).'));
    
    pl_01cell = [ 0 1 pl_cell(:).' ];
    names = [ '0' '1' cellfun(@(pl) pl.blk.names(1), pl_cell) ];
    
    [im,indices] = ismember(paffmat_names, names);
    
    if ~all(im)
        missing = strjoin(names(~im),', ');
        error('LFR variables %s must appear in the second argument of PAffineMatrix.lfr(~,pl_cell)', missing);
    end
    
    Xi_cell = cellfun(@(ind) { N.Im * pl_01cell{ind} }, num2cell(indices));   
    
else
    
    Xi_cell = cellfun(@(c) { N.Im * channel2lfr_simple(c) }, num2cell(N.channels));

end
    
Xi = vertcat(Xi_cell{:});

syslfr = N.Theta_left * Xi;

end

function c = channel2lfr_simple(c)
    if isempty(symvar(c))
        c = lfr(double(c));
    elseif isscalar(symvar(c))
        c = lfr(char(c));
    end
end

function test1
%%
% 2020.03.31. (március 31, kedd), 00:35

syms p1 x2 real

[~,p_lfr_cell] = pcz_generateLFRStateVector('p',[-2 3 ; -3 5],[-2 2 ; -3 3]);
[~,x_lfr_cell] = pcz_generateLFRStateVector('x',[-12 13 ; -31 15],[-12 12 ; -13 13]);

paffmat = PAffineMatrix([ 1 2 3 4 5 6 ], [1;p1;x2])

G = plfr(paffmat,[ p_lfr_cell x_lfr_cell ])

OK_IF_ZERO = sym(paffmat) - sym(G)

end

function test2
%%
% 2020.03.31. (március 31, kedd), 00:35

He = he;

np = 2;
p_lim = repmat([4,5],[np,1]);
dp_lim = repmat([-4,4],[np,1]);

P_generate_symvars(nx,np,nu,ny);
[p_lfr,p_lfr_cell] = pcz_generateLFRStateVector('p',p_lim,dp_lim);
[dp_lfr,dp_lfr_cell] = pcz_generateLFRStateVector('dp',p_lim);
pdp_lfr_cell = [ p_lfr_cell dp_lfr_cell ];

m = 2;
N = PAffineMatrix(cellfun(@(i){ round(He(0.49*rand(m))) }, num2cell(0:2*np)), [1 ; pdp])

N_lfr = plfr(N,pdp_lfr_cell)

OK_IF_ZERO = sym(N) - sym(N_lfr)

end

function test3
%%

[p,pl] = pcz_generateLFRStateVector('p',[2 3 ; 3 5],[-2 2 ; -3 3]);

PI = [
    eye(2)
    eye(2)*pl{1}
    eye(2)*pl{2}
    eye(2)*pl{1}^2
    eye(2)*pl{1}*pl{2}
    eye(2)*pl{2}^2
    ];

[N,~,~,N_err] = P_affine_annihilator_for_LFR(PI,p,'lims',[-1 1 ; -1 1],'sym',1);
pcz_lfrzero_report(minlfr(N*PI), N_err, 'Annihilator N');

sym(N)


end