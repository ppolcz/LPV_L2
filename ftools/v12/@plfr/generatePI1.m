function PI1 = generatePI1(pLFR)
%%

s1 = pLFR.ny;
s2 = pLFR.nu;
m = s2+pLFR.m1;

[~,~,M11,M12] = pcz_split_matrix(eye(m),[s2 pLFR.m1],[s2 pLFR.m1]);

PI1 = plfr(M11,M12,pLFR.C,pLFR.D,pLFR.blk);
PI1 = PI1.set_vars(pLFR.subsvars(:));

end
