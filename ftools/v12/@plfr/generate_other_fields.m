function pLFR = generate_other_fields(pLFR)
%%

% TODO: nem lenne-e jobb dependent var-ket kezelni?
pLFR.M = [pLFR.A,pLFR.B;pLFR.C,pLFR.D];
pLFR.I = eye(pLFR.m1);        

end
