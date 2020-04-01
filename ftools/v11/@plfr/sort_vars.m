function [pLFR,new_vars_] = sort_vars(pLFR)
%%
% Created:  2020.03.29. (március 29, vasárnap), 12:19

names = pLFR.names(~strcmp('1',pLFR.names));

if isempty(names)
    pLFR = pLFR.set_vars([]);
    return
end

% E.g.
% parsed{3}.name = 'p'
% parsed{3}.nr = 1
%
% { [name: 'dp', nr: 1] [name: 'dp', nr: 2] [name: 'p', nr: 1] ... }
parsed = cellfun(pLFR.Parse_VarName, names, 'UniformOutput', false);

% { 'dp' 'dp' 'p' 'p' 'x' 'a' }
varnames = cellfun(@(p) {p.name}, parsed);

% [ 20 20 10 10 1 Inf ]
order_nrs = cellfun(@(vn) pLFR.ordernr__(vn), varnames);

% [ 6 5 3 4 1 2 ]
[~,ind] = sort(order_nrs);

% [ x1 p1 p2 dp1 dp2 a ]
new_vars_ = pLFR.symvars(ind);

pLFR = pLFR.set_vars(new_vars_);

end
