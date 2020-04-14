function method0_grid_ltiwc_Hinf(modelname,A_fh,B_fh,C_fh,D_fh,p_lim,varargin)
%% method0_grid_ltiwc_Hinf
%  
%  File: method0_grid_ltiwc_Hinf.m
%  Directory: 8_published/LPV_L2/workspace
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2020. April 14. (2019b)
%

%%

args.res = 10;
args = parsepropval(args,varargin{:});

np = size(p_lim,1);

resolution = args.res;
if isscalar(resolution)
    resolution = repmat(resolution,[np 1]);
end

[A_fh,B_fh,C_fh,D_fh] = helper_fh2fh(A_fh,B_fh,C_fh,D_fh);


% Generate grid
lspace = cellfun(@(o) {linspace(o{:})}, num2cell(num2cell([p_lim resolution]),2));
pp = cell(1,np);
[pp{:}] = ndgrid(lspace{:});
pp = cellfun(@(a) {a(:)'}, pp);
pp = num2cell(num2cell(vertcat(pp{:})),1);

Hinf = @(p) norm(ss(A_fh(p{:}),B_fh(p{:}),C_fh(p{:}),D_fh(p{:})),Inf);

Hinf_all = cellfun(@(p) Hinf(p), pp);

[wcg,wci] = max(Hinf_all);

wcp_str = cellfun(@(n,i) {sprintf('p%d = %s',i,num2str(n))}, pp{wci}, num2cell(1:np)');

pcz_dispFunction2('Worst case LTI Hinf norm: %g', wcg)
pcz_dispFunction2('Worst case parameter values: \n  %s', strjoin(wcp_str,[newline '  ']))

store_results('Results_All.csv', modelname, wcg, 0, 0, 0, ...
    strjoin(wcp_str,', '), 'Grid-based worst case LTI Hinf analysis (frozen parameters)')


end