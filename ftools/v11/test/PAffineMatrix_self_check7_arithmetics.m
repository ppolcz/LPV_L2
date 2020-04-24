%% PAffineMatrix_self_check7_arithmetics
%  
%  File: PAffineMatrix_self_check7_arithmetics.m
%  Directory: 7_ftools/ftools/v11/test
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2019. November 05. (2019a)
%

%% minus

syms p q r real
vars = [1;p;q;r];

X = PAffineMatrix(round(0.5*randn(4,12)),vars);
Y = PAffineMatrix(round(0.5*randn(4,12)),vars);

X.Sym = 1;
Y.Sym = 1;

Z = X - Y;

Z.Sym = 1;

pcz_symzero(sym(X) - Y - Z, 'X - Y = Z')

%% mtimes

Theta = [
    1 0 0 0 0 1
    0 1 0 1 0 0
    ];

syms p q real

A = PAffineMatrix(Theta, [1;p;q]);
A.Sym = 1

A_sym = A.Sym;

A = A * [1 ; 1]
A = [ 2 1 ] * A

display(A.Theta, 'A.Theta')

pcz_symzero(A.Sym - [2 1] * A_sym * [1 ; 1], 'left and right matrix multiplication')


%% plus1

syms p q r real
vars = [1;p;q;r];

X = PAffineMatrix(round(2*randn(4,12)),vars);
Y = PAffineMatrix(round(2*randn(4,12)),vars);

X.Sym = 1
Y.Sym = 1

Z = X + Y;

Z.Sym = 1

pcz_symzero(sym(X) + Y - Z, 'X - Y = Z')

%% plus2

syms p real
tmp = PAffineMatrix(ones(3,3),p)
tmp + eye(3)
tmp = tmp + eye(3)
tmp.Sym = 1
