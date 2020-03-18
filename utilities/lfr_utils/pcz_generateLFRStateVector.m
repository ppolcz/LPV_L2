function [x,x_cell] = pcz_generateLFRStateVector(name,dim_bounds,varargin)
%% pcz_generateLFRStateVector
%  
%  File: pcz_generateLFRStateVector.m
%  Directory: 2_demonstrations/lib/matlab
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2019. April 06.
%

if nargout > 0
    props.GenerateVars = 0;
else
    props.GenerateVars = 1;
end
props = parsepropval(props,varargin{:});

%%

if isscalar(dim_bounds)
    dim = dim_bounds;
    bounds = [ -ones(dim,1) ones(dim,1) ];
else
    bounds = dim_bounds;
    dim = size(bounds,1);
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

x_names = cellfun(@(i) {[name num2str(i)]}, num2cell(1:dim));

x_cell = cellfun(@(i,name_i) { lfr(name_i,'ltisr',1,bounds(i,:),'minmax') }, num2cell(1:dim), x_names);
x = vertcat(x_cell{:});

if props.GenerateVars == 1
    for i = 1:dim    
        assignin('caller', x_names{i}, x_cell{i});
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