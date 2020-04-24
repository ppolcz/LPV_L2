%% 
%  File: zygote_self_check4_quadratic_with_affine_then_subspace.m
%  Directory: 1_PhD_projects/00_my_toolboxes/FinslerTools/v11/demos
%  Author: Peter Polcz (ppolcz@gmail.com) 
% 
%  Created on 2018. October 23.
% 
%%

global SCOPE_DEPTH VERBOSE LATEX_EQNR 
SCOPE_DEPTH = 0;
VERBOSE = 1;
LATEX_EQNR = 0;

TMP_sACiDKwRMsxuOMwKXOKf = pcz_dispFunctionName;

n = 2;
d_n = 1;
P_generate_symvars_v5(n,d_n);

Pi = monomials(x,2);
p = numel(Pi);

Theta = cellfun(@(i){ sdpvar(p,p,'symmetric') }, num2cell(0:d_n));
P = PAffineMatrix(Theta,[1;d]);
P = P.to_symTheta('P');
P.Sym = 1


V = x1^4 + d1*x2^4;


V_all = expand(Pi.' * P * Pi)

c = coeffs(V_all - V, xd);

p = symvar(P.Theta);
Theta_fh = matlabFunction(P.Theta, 'vars', {transpose(p)})

[A,B] = equationsToMatrix(c == 0, p);

A = double(A);
B = double(B);


p0 = round(A\B,4);

Ker_A = null(A);
dim_Ker = size(Ker_A,2);


if dim_Ker == 0
    P_span = p0;
else
    P_span = Ker_A + p0(:,ones(1,dim_Ker));
end
P_span = round(P_span,10);

dim_P = size(P_span,2);

alpha = sdpvar(dim_P,1);
Theta_P_ss_cell = cellfun(@(p) {Theta_fh(p)}, num2cell(P_span,1));
Theta_P_ss = PAffineMatrix(Theta_P_ss_cell,alpha);
Theta_P_ss_sym = Theta_P_ss.to_symTheta('a');

P = PAffineMatrix(Theta_P_ss.Matrix, [1;d]);
P = P.set_subsvars(xd);

P_sym = PAffineMatrix(Theta_P_ss_sym.Matrix, [1;d]);
P_sym.Sym = 1;

V_param = simplify(Pi' * P_sym.Sym * Pi);

N = P_affine_annihilator(Pi,x,'sym',1);
N = N.set_subsvars(xd);

L = sdpvar(N.m,N.q,'full');

xd_lim = [
    -1 1
    -1 1
    1 1.2
    ];

XD_v = P_ndnorms_of_X(xd_lim);

CONS = sum(alpha) == 1;
% CONS = [];
for i = 1:size(XD_v,1)
    xd_num = XD_v(i,:)';
    
    N_num = N(xd_num);
    
    CONS = [ CONS
        P(xd_num) + L*N_num + N_num'*L' - 0.001*eye(size(P)) >= 0 
        ];
end

spdopts = sdpsettings('solver','sedumi');

sol = optimize(CONS,[],spdopts)

P = value(P);

simplify(Pi' * P.Matrix * Pi / sum(value(alpha)))

pcz_dispFunctionEnd(TMP_sACiDKwRMsxuOMwKXOKf);
