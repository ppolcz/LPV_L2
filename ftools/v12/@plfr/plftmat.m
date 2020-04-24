function G_plftmat = plftmat(G)
%% plftmat
%  
%  File: plftmat.m
%  Directory: 7_ftools/ftools/v11/@plfr
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2020. March 30. (2019b)
%

%%

deltai_fh_cell = cellfun(@(s) { matlabFunction(s,'vars',G.symvars) }, num2cell(diag(G.Delta)));

n = cell(G.varnames).';
b = num2cell(G.bounds,2);
r = num2cell(G.rbounds,2);

tvr_cell = cellfun(@(n,b,r) { tvreal(n,b,r) }, n, b, r);

deltai_tvr_cell = cellfun(@(d) { d(tvr_cell{:}) }, deltai_fh_cell);

Delta_plftmat = blkdiag(deltai_tvr_cell{:});

G_plftmat = G.A + G.B * ( (G.I - Delta_plftmat * G.D) \ Delta_plftmat ) * G.C;

end


function test1
%%

syms a1 d1 p1 p2 p3 t1 t2 real

vars = [a1 d1 p1 p2 p3 t1 t2].';

p_lim = [-1 1 ; 0 1 ; 1 2];
dp_lim = [-1 1; -2 2 ; -6 6];

d_lim = [9 10];
Dd_lim = [-1 1];

a_lim = [1 10];
da_lim = [-10 10];

[p_lfr, p_lfr_cell] = pcz_generateLFRStateVector('p',p_lim,dp_lim);
[d_lfr, d_lfr_cell] = pcz_generateLFRStateVector('d',d_lim,Dd_lim);
[a_lfr, a_lfr_cell] = pcz_generateLFRStateVector('a',a_lim,da_lim);

H = plfr([p_lfr;d_lfr;a_lfr]);

vars_bounds_rbounds = [ H.symvars.' H.bounds H.rbounds ]

ELVART = [
    a1 a_lim da_lim
    d1 d_lim Dd_lim
    [p1 p2 p3].' p_lim dp_lim
    ];

ZERO = vars_bounds_rbounds - ELVART;



A_fh = @(p1,p2,p3,a1,d1) [
    -3+p1       3+p1/(p3^2 + 0.5*d1 + 1)     0.1*p3   p3^2*p2
    0          -1-p2^2                       5        0
    -1/(5-p2)  0                             -4+p1    1/a1
    0          0.1                           0        -5+1/(p1 + 2)
    ];
    

G = plfr(A_fh(p_lfr_cell{:},a_lfr_cell{:},d_lfr_cell{:}));


G = G.set_vars(vars)


G_plftmat = plftmat(G)

end