function [p,p_cell] = pcz_generateGSStateVector(name,dim_or_lim,dp_lim,varargin)
%% pcz_generateLFRStateVector
%
%  File: pcz_generateLFRStateVector.m
%  Directory: 2_demonstrations/lib/matlab
%  Author: Peter Polcz (ppolcz@gmail.com)
%
%  Created on 2019. April 06.
%

%%

if isscalar(dim_or_lim)
    dim = dim_or_lim;
    p_lim = [ -ones(dim,1) ones(dim,1) ];
    dp_lim = [ -ones(dim,1) ones(dim,1) ];
else
    p_lim = dim_or_lim;
    dim = size(p_lim,1);

    if nargin < 3 || isempty(dp_lim)
        dp_lim = [ -ones(dim,1) ones(dim,1) ]*Inf;
    end
end

p_nom = p_lim * [1 ; 1]/2;


p_nom_c = num2cell(p_nom,2);
p_lim_c = num2cell(p_lim,2);
dp_lim_c = num2cell(dp_lim,2);

names = cellfun(@(i) {[name num2str(i)]}, num2cell(1:dim));

p_cell = cellfun(@gss,names(:),p_nom_c,p_lim_c,dp_lim_c,'UniformOutput',false);
p = vertcat(p_cell{:});

end

function test1
%%

[p,p_cell] = pcz_generateGSStateVector('x',[-1 1 ; -2 2 ; -3 3])

[p,p_cell] = pcz_generateGSStateVector('x',4)


end