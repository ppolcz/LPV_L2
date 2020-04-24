function pLFR = generate_LFR(pLFR,blk)
%%

bounds = [];
if size(blk,1) == 2
    bounds = blk;
    blk = [];
end
    

% Elotte_Gx = pLFR.A + pLFR.B/(I - pLFR.Delta*pLFR.D)*pLFR.Delta*pLFR.C;

[delta,sigma] = sort(diag(pLFR.Delta));

pLFR.Delta = diag(delta);

Is = pcz_permat(sigma);

pLFR.A = pLFR.A;
pLFR.B = pLFR.B*Is;
pLFR.C = Is'*pLFR.C;
pLFR.D = Is'*pLFR.D*Is;
pLFR.Delta = Is' * pLFR.Delta * Is;

% Utana_Gx = A + B/(I - Delta*D)*Delta*C;
% pcz_symzero_report(tGx - Gx, 'new G(x) = G(x) before permutation');

[p,cumnr] = unique(delta);

dim = diff([ cumnr ; size(pLFR.Delta,1)+1 ])';

% dim
% pLFR.bounds'

lfrtbx_blk.names = cellfun(@(pi) {char(pi)}, num2cell(p.'));

if ~isempty(blk)

    % We have already a desc matrix, but we need to update it.

    % lfrtbx_blk.names
    % blk.names

    idx = zeros(1,numel(lfrtbx_blk.names));
    for i = 1:numel(lfrtbx_blk.names)
        idx(i) = find(strcmp(blk.names,lfrtbx_blk.names{i}));
    end

    lfrtbx_blk.desc = blk.desc(:,idx);

    lfrtbx_blk.desc(1,:) = dim;
    lfrtbx_blk.desc(2,:) = dim;

else

    np = numel(symvar(pLFR.Delta));
    
    if size(bounds,2) == np
        bounds = [
            bounds
            [0.5 0.5]*bounds
            [-Inf;Inf]*ones(1,np)
            ];
    elseif size(bounds,2) == 2*np
        bounds = [
            bounds(:,1:np)
            [0.5 0.5]*bounds(:,1:np)
            bounds(:,np+1:2*np)
            ];
    end
    
    % If the bounds does not contain the first 0 0 values for the
    % constant part of the Delta block
    if p(1) == 1 && np == numel(p)-1
        bounds = [ zeros(5,1) bounds ];
    end

    % Delta interface for LFR Toolbox
    lfrtbx_blk.desc = [
        % row dimension
        dim
        % column dimension
        dim
        % real(1)/complex(0)
        dim*0 + 1
        % scalar(1)/full(0)
        dim*0 + 1
        % linear(1)/nonlinear(0)
        dim*0 + 1
        % time-invariant(1)(memoryless if nonlinear)/time-varying(0)
        dim*0 + 1
        % bound informations (given as intervals: 1,2)
        dim*0 + 1
        dim*0 + 2
        % b.i.: intervals
        bounds
        ];
end

pLFR.lfrtbx_obj = lfr(pLFR.D,pLFR.C,pLFR.B,pLFR.A,lfrtbx_blk);

end