%% plfr_test
%
%  File: plfr_test.m
%  Directory: 7_ftools/ftools/v11/test
%  Author: Peter Polcz (ppolcz@gmail.com)
%
%  Created on 2019. November 06. (2019a)
%


G_reset

%% test1

lfrs x b dp [1,2,3] [7,8,9] [4.2,5.1,6.5]

syslfr = [ x^2+b ; b+1 ; x/(dp^2 + 1) ];

pLFR = plfr(syslfr)

pcz_symzero_report(pLFR(1,2,3) - [1+3 ; 3+1 ; 1/(2^2 + 1)])

%% TEST compatibility with P_LFR_reduction_numeric

[x,x_cell] = pcz_generateLFRStateVector('x',4);

Sf = 10;
Kr = 1;
Y = 0.5;

f_fh = @(x1,x2,x3,F) [
    Kr*x1*x2 - x1/x3*F
    -Kr/Y*x1*x2 + (Sf-x2)/x3*F
    F
    0
    ];

pLFR = plfr(f_fh(x_cell{:}));

Delta = pLFR.Delta
lfr_desc = pLFR.desc

pLFR_min = P_LFR_reduction_numeric(pLFR)

sym(pLFR)

A = [pLFR.A pLFR.B]
PI = sym(pLFR.generatePI)

sym(pLFR_min)

A_min = [pLFR_min.A pLFR_min.B]
PI_min = sym(pLFR_min.generatePI)

pcz_symzero_report(A*PI - A_min*PI_min, 'Reduced LFR gives the same?')

pcz_symzero_report(A*PI - sym(plfr(pLFR_min.lfrtbx_obj)), 'Reduced LFR u.a. (meg ha az lfrtbx_obj-et veszem, akkor is!')

%% test2
lfrs a b c [1,2,3] [7,8,9] [4.2,5.1,6.5]

syslfr = [ a^2+b ; b+1 ; a/(c^2 + 1) ];

[D,C,B,A] = lfrdata(syslfr);

pLFR = plfr([A B ; C D],syslfr.blk)

pLFR(1,2,3)

a = 1;
b = 2;
c = 3;

ELVART = [ a^2+b ; b+1 ; a/(c^2 + 1) ]

%% test3
pcz_generateSymStateVector(4,'x','real');

s1 = 2;
s2 = 3;
m1 = 9;

M = round(0.5*randn(s1+m1,s2+m1));
Delta = diag([x4 x2 x3 x2 x4 x4 x1 x3 x3]);

p_lims = [
    -1 1
    -1 1
    -1 2
    -4 4
    ];

pLFR = plfr(M,Delta,p_lims) %#ok<NOPRT>

pLFR(1,2,3,5)

%% test6
pcz_generateSymStateVector(4,'x','real');

s1 = 2;
s2 = 3;
m1 = 9;

[A,B,C,D] = pcz_split_matrix(round(0.5*randn(s1+m1,s2+m1)),[s1 m1],[s2 m1]);
Delta = diag([x4 x2 x3 x2 x4 x4 x1 x3 x3]);

p_lims = [
    -1 1
    -1 1
    -1 2
    -4 4
    ];

pLFR = plfr(A,B,C,D,Delta,p_lims) %#ok<NOPRT>

pLFR([1,2,3,4]')
pLFR(1,2,3,4)

%% test6_const
pcz_generateSymStateVector(4,'x','real');

s1 = 2;
s2 = 3;
m1 = 9;

[A,B,C,D] = pcz_split_matrix(round(0.5*randn(s1+m1,s2+m1)),[s1 m1],[s2 m1]);
Delta = diag([1 x2 x3 x2 1 x4 x1 x3 x3]);

p_lims = [
    -1 1
    -1 1
    -1 2
    -4 4
    ];

pLFR = plfr(A,B,C,D,Delta,p_lims) %#ok<NOPRT>

pLFR([1,2,3,4]')
pLFR(1,2,3,4)

%% test6_no_bounds
pcz_generateSymStateVector(4,'x','real');

s1 = 2;
s2 = 3;
m1 = 9;

[A,B,C,D] = pcz_split_matrix(round(0.5*randn(s1+m1,s2+m1)),[s1 m1],[s2 m1]);
Delta = diag([x4 x2 x3 x2 x4 x4 x1 x3 x3]);

pLFR = plfr(A,B,C,D,Delta,[]) %#ok<NOPRT>

pLFR([1,2,3,4]')
pLFR(1,2,3,4)

%% test_EMPTY -- ERRE NEM MUKODIK, DE NEM IS KELL

% pLFR = plfr;


%% Constructor test - van 1-es block is

[x,x_cell] = pcz_generateLFRStateVector('x',4);

Sf = 10;
Kr = 1;
Y = 0.5;

f_fh = @(x1,x2,x3,F) [
    Kr*x1*x2 - x1/x3*F
    -Kr/Y*x1*x2 + (Sf-x2)/x3*F
    F
    0
    ];

N_PI_fh = @(I,x1,x2,x3,F) [
    I
    I*x1
    I*x2
    I*x3
    I*F
    I*(x1/Y - x2 + Sf)^(-1)
    ];
% -------------

n = numel(x_cell);

syslfr = f_fh(x_cell{:});

pLFR = P_lfrdata(syslfr);

% case 1
pLFR1 = plfr(syslfr)

% case 2
pLFR2 = plfr(pLFR1.M,pLFR1.blk)

% case 3
pLFR3 = plfr(pLFR1.M,pLFR1.Delta,pLFR1.bounds)

% case 5
pLFR5 = plfr(pLFR1.A,pLFR1.B,pLFR1.C,pLFR1.D,pLFR1.blk)

% case 6
pLFR6 = plfr(pLFR1.A,pLFR1.B,pLFR1.C,pLFR1.D,pLFR1.Delta,pLFR1.bounds)

% case 7
pLFR7 = plfr(pLFR1.A,pLFR1.B,pLFR1.C,pLFR1.D,pLFR1.Delta,[],pLFR1.blk)

pcz_symzero(sym(pLFR1) - sym(pLFR2), 'LFR1 == LFR2')
pcz_symzero(sym(pLFR1) - sym(pLFR3), 'LFR1 == LFR3')
% pcz_symzero(sym(pLFR1) - sym(pLFR4), 'LFR1 == LFR4')
pcz_symzero(sym(pLFR1) - sym(pLFR5), 'LFR1 == LFR5')
pcz_symzero(sym(pLFR1) - sym(pLFR6), 'LFR1 == LFR6')
pcz_symzero(sym(pLFR1) - sym(pLFR7), 'LFR1 == LFR7')


%% Reciprok


lfrs p
G = 1/p;

[A,B,C,D,Delta] = P_lfrdata_v1(G);

pLFR = plfr(G)

%% Order of variables
% 2020.03.29. (március 29, vasárnap), 12:14

lfrs x1 x2 p1 p2 dp1 a b c d x3
pLFR = plfr(x1+x2+p1+p2+dp1+a+b+c+d+x3)


%% Diff

[p,p_cell] = pcz_generateLFRStateVector('p',4);
[dp,dp_cell] = pcz_generateLFRStateVector('dp',4);

pLFR = plfr([ 
    p(1) + p(2)^2
    p(3)*p(4)
    ]);

diff(pLFR,p_cell,dp_cell)

%%

P_init(11)
G_reset

x = pcz_generateGSStateVector('x',1);
p = pcz_generateGSStateVector('p',1);

xl = pcz_generateLFRStateVector('x',1);
pl = pcz_generateLFRStateVector('p',1);

fh = @(x,p) [ 
    p^2      1/(x^2+p+5)     1/(p^2 + p + 1) 
    1/(x+10) x^3/(x^4+p^2+1) x^5*p^7/(x+p^2+10) ];

pp = fh(x,p);

ppl = plfr(fh(xl,pl));

N = 20*1000;

TMP_UIWKgQgOJnYfrgPwpoUB = pcz_dispFunctionName('gss object');
eval(pp, {'x1' 'p1'}, num2cell(rand(N,2)));
pcz_dispFunctionEnd(TMP_UIWKgQgOJnYfrgPwpoUB);

TMP_UIWKgQgOJnYfrgPwpoUB = pcz_dispFunctionName('plfr object');
ppl.val(rand(N,2)); 
pcz_dispFunctionEnd(TMP_UIWKgQgOJnYfrgPwpoUB);

%%

p_lim = [
    -1 2
    -1 2
    0 2
    ];

dp_lim = [
    -10 10
    -1 1
    -5 5
    ];

[p_gss,p_gss_cell] = pcz_generateGSStateVector('p',p_lim,dp_lim);

A_fh = @(p1,p2,p3) [
    -3+p1       3+p1/(p3^2 + 0.5*p1 + 1)     0.1*p3   p3^2*p2
    0          -1-p2^2                       5        0
    -1/(5-p2)  0                             -4+p1    0
    0          0.1                           0        -5+1/(p1 + 2)
    ];

A_gss = A_fh(p_gss_cell{:});

[A_samples,samples] = dbsample(A_gss,10)
