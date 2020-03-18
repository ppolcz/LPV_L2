function [N] = P_affine_annihilator(Pi, w, x, varargin)
%% P_annihilator_linsolve_v9
% 
%  File: P_affine_annihilator.m
%  Directory: 1_PhD_projects/00_my_toolboxes/FinslerTools/v11
%  Author: Peter Polcz (ppolcz@gmail.com) 
% 
%  Created on 2018. October 22.

%%


TMP_qftrrQBBOgGsrDSkIUkp = pcz_dispFunctionName;


%%

% Compatibility with the older version, where the second vector of
% variables are not given:
% 2018.09.27. (szeptember 27, csütörtök), 12:24
narginchk(2,15);
if nargin < 3 || isempty(x)
    x = w;
elseif ischar(x)
    varargin = [ x varargin];
    x = w;
end

args.caption = '';
args.symbolic = 0;
args.whitening = 0;
args.precision = 10;
args.returncell = 0;
args.tol = 1e-10;
args = parsepropval(args,varargin{:});

%% Annihilator generalas

% w = transpose(symvar(Pi));

s = numel(w);
n = numel(x);
m = numel(Pi);

b = sym('B__', [1,m]);
a = sym('A__', [n,m]);

% Egy teljesen altalanos affin annihilator sor
% 2018.09.27. (szeptember 27, csütörtök), 12:25
row_sym = x.' * a + b;

% A szamlalo nulla kell legyen
[num,den] = numden(row_sym*Pi);

% Numerikus szamitasi hibak csokkentese erdekeben a szamlalot elosztom
% a nevezo legnagyobb egyutthatojaval
c = max(double(coeffs(den)));
num = collect(num / c(1), w);

vars = [b ; a];
% num
% symvar(num)
% w
[A,B] = equationsToMatrix(coeffs(num,w) == 0, vars(:));

% 2018.09.27. (szeptember 27, csütörtök), 12:25
assert(isempty(symvar(A)), ...
    '[A,] = equationsToMatrix, A should not contain symbolic variables, symvar(A) = %s', ...
    char(symvar(A)))

assert(isempty(symvar(B)), ...
    '[,B] = equationsToMatrix, B should not contain symbolic variables, symvar(B) = %s', ...
    char(symvar(B)))

A = double(A);
assert(all(double(B) == 0), ...
    'Matrix B in the linsolve should be zero!');

N1_coeffs = pcz_vecalg_null(A,args.tol)';

N = PAffineMatrix(N1_coeffs, [1;x], w);
N.caption = args.caption;



% POST PROCESSING



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [1] Numberical Quick check I.

Pi_fh = matlabFunction(Pi, 'vars', {w});
pcz_fhzero_report(@(w) N(w) * Pi_fh(w), N.subsvars,...
    [ 'Numerical check. ' N.caption ])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [2] beautify annihilator
if args.symbolic || args.whitening
    
    TMP_iHMHxCpMjJjtGOeWOzwx = pcz_dispFunctionName('beautify annihilator + !roundings!');

    % itt tortenik a beautifyingolas
    N2_coeffs = rref(N1_coeffs);
    
    % itt tortenik a roundingolas !!
    % N3_coeffs = round(N2_coeffs, args.precision + 2);
    N3_coeffs = double(pcz_find_recdec(N2_coeffs, 'tolerance', 10^(-args.precision)));
    
    N.Theta = N3_coeffs;
    N = N.generate_symbolic;    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % [2.1] Numerical Quick check II.
    
    pcz_fhzero_report(@(w) N(w) * Pi_fh(w), N.subsvars,...
        [ 'Numerical check. ' N.caption ])

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % [2.2] Symbolical check
    
    pcz_symzero_report(simplify(N.Sym * Pi),...
        [ 'Symbolical check. ' N.caption ])

    pcz_dispFunctionEnd(TMP_iHMHxCpMjJjtGOeWOzwx);
    clear TMP_iHMHxCpMjJjtGOeWOzwx    

end

%%
pcz_dispFunctionEnd(TMP_qftrrQBBOgGsrDSkIUkp); 

end


%% Self check

%#ok<*DEFNU>

function self_check1
%%
    global SCOPE_DEPTH VERBOSE 
    SCOPE_DEPTH = 0;
    VERBOSE = 1;

    P_generate_symvars_v5(4,2);
    
    Pi = [
        x
        x*d1
        x*d2
        x1*x2*d2
        x1*x2^2*d2
        ];
    
    w = xd;
    x = d;
    
    varargin = {'sym', 1, 'caption', 'Nb: Self check.' };
    
    N = P_affine_annihilator(Pi, xd, d, varargin{:});
    
    N.sym
        
end

