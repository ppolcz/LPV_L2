function [x,x_cell] = pcz_generateLFRStateVector(name,dim_or_lim,dp_lim,varargin)
%% pcz_generateLFRStateVector
%
%  File: pcz_generateLFRStateVector.m
%  Directory: 2_demonstrations/lib/matlab
%  Author: Peter Polcz (ppolcz@gmail.com)
%
%  Created on 2019. April 06.
%

props.GenerateVars = 0;
props = parsepropval(props,varargin{:});

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

% % Delta interface for LFR Toolbox
% blk.desc = [
%     % 1: row dimension
%     % 2: column dimension
%     % 3: real(1)/complex(0)
%     % 4: scalar(1)/full(0)
%     % 5: linear(1)/nonlinear(0)
%     % 6: time-invariant(1)(memoryless if nonlinear)/time-varying(0)
%     % 7,8: bound informations (given as intervals: 1,2)
%     ones(7,dim)    % 1-7
%     ones(1,dim)*2  % 8
%     % b.i.: intervals
%     bounds'
%     % b.i.: nominal values
%     [0.5 0.5]*bounds'
%     ];
% blk.names = cellfun(@(i) { [ name num2str(i) ] }, num2cell(1:dim));
%
% A = zeros(dim,1);
% B = eye(dim);
% C = ones(dim,1);
% D = zeros(dim);
%
% x = lfr(D,C,B,A,blk);

desc = [
    %  1: row-dimensions of blocks
    %  2: column-dimensions of blocks
    %  3: real(1) / complex(0) block types
    %  4: scalar(1) / full(0) block types
    %  5: linear(1) / nonlinear (0) block types
    %  6: time-inv.(1) / time-var.(0) block types
    %  7: min/max(1) / sector(2) / freq. dependent(>2) bounds, vagy 0
    1*ones(7,dim)

    % 8: min/max(2) / sector(1) / freq. dependent(>2) bounds, vagy 0
    2*ones(1,dim)

    %  9: minumum values of bounds
    % 10: maximum values of bounds
    p_lim'

    % 11: nominal value
    [0.5 0.5]*p_lim'
    
    % 2020.03.29. (március 29, vasárnap), 20:35
    % -------------------------------------
    % 12: minimum rate bound
    % 13: maximum rate bound
    dp_lim'
    ];
    
desc_cell = num2cell(desc,1);
names_cell = cellfun(@(i) {[name num2str(i)]}, num2cell(1:dim));

blk_cell = cellfun( @(d,n) {struct('names',{{n}},'desc',d)}, desc_cell, names_cell);

x_cell = cellfun( @(b) {lfr(0,1,1,0,b)}, blk_cell );

%x_cell = cellfun(@(i,name_i) { lfr(name_i,'ltisr',1,p_lim(i,:),'minmax') }, num2cell(1:dim), names_cell);

x = vertcat(x_cell{:});

if props.GenerateVars == 1
    for i = 1:dim
        assignin('caller', names_cell{i}, x_cell{i});
    end

    assignin('caller', name, x);
    assignin('caller', [name '_cell'], x_cell);
end

end

function test1
%%

[p,p_cell] = pcz_generateLFRStateVector('x',[-1 1 ; -2 2 ; -3 3])

[p,p_cell] = pcz_generateLFRStateVector('x',4)


end
