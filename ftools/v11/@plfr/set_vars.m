function pLFR = set_vars(pLFR,vars_)
%%

nem_szerepel = setdiff(pLFR.symvars,vars_);

if ~isempty(nem_szerepel)
    error('plfr: Variables `%s` must be included!', ...
        strjoin(cellfun(@(x) {char(x)}, num2cell(nem_szerepel)),', '))
end

pLFR.subsvars = vars_(:);

delta = diag(pLFR.Delta);

if isempty(delta)
    pLFR.delta_fh = @(xp) zeros(0,size(xp,2));
else
    ZERO = sym('ZERO');
    delta(delta == 1) = 1 + ZERO;

    pLFR.delta_fh = matlabFunction(delta,'vars',[ pLFR.subsvars(:) ; ZERO ]);
end

% pLFR.deltai_fh_cell = cellfun(@(s) { matlabFunction(s,'vars',G.symvars) }, num2cell(diag(G.Delta)));

end