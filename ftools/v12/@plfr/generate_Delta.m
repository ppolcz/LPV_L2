function pLFR = generate_Delta(pLFR)
%%
% Requires a valid lfrtbx_obj


% 2020.03.30. (március 30, hétfő), 03:08
if size(pLFR.desc,1) <= 11
    pLFR.lfrtbx_obj.blk.desc(12,pLFR.desc(10,:) ~= 0) = -Inf;
    pLFR.lfrtbx_obj.blk.desc(13,pLFR.desc(10,:) ~= 0) = Inf;
end

s = numel(pLFR.names);
c = cell(s, 1);
for k = 1:s
    c{k} = sym(pLFR.names{k}) * eye(pLFR.desc(1:2,k)');
end
pLFR.Delta = blkdiag(c{:});


end
