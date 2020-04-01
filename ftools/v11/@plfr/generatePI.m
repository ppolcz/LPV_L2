function PI = generatePI(pLFR)
%%

s1 = pLFR.ny;
s2 = pLFR.nu;
m = s2+pLFR.m1;

[M11,M12] = pcz_split_matrix(eye(m),m,[s2 pLFR.m1]);

PI = plfr(M11,M12,pLFR.C,pLFR.D,pLFR.blk);
PI = PI.set_vars(pLFR.subsvars(:));

end