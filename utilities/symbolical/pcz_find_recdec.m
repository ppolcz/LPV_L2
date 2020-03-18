function [A_sym,A_int,A_num,A_den,A_diff,A_not_Good] = pcz_find_recdec(A, varargin)
%% 
%  File: pcz_find_recdec.m
%  Directory: utilities/symbolical
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2019. June 07. (2019a)
%  Modified on 2019. November 20. (2019a)
%

%%
 
args.tolerance = 1e-10;
args.symbolic = 1;

% 2019.11.20. (november 20, szerda), 14:21
args.maxden = 1000;
args.maxit = 6;
args.step = 10;

% 2019.11.20. (november 20, szerda), 14:30
args.structout = false;

%%

args = parsepropval(args,varargin{:});

%%

A_ = A;

r = @(A) reshape(A,size(A_));

A = A(:);

ret1 = A;
ret2 = 1;

A_int = floor(A);
A_num = A - A_int;
A_den = ones(size(A));
A_diff = zeros(size(A));

good = A_num < args.tolerance;
A_num(good) = 0;

if prod(good) > 0
    % [ret,varargout{:}] = prepare_return_values(A_integer_part,A,1,args);

    max_denominator = max(A_den);

    A_sym = r(sym(A));
    A_int = r(A_int);
    A_num = r(A_num);
    A_den = r(A_den);
    A_diff = r(A_diff);
    
    A_int = prepret;
    return
end

% 2019.11.20. (november 20, szerda), 14:21
min_denominator = 2;
max_denominator = args.maxden;
for i = 1:args.maxit
    
    denominators = min_denominator:max_denominator;

    A_not_Good = A_num(~good);
    
    s = numel(A_not_Good);
    p = numel(denominators);

    % Find potential numerators, rounded to the closest integer:
    numerators = round(A_not_Good * denominators);
    
    % We need to get rid of zeros.
    % zeroIndexes = numerators == 0;
    % numerators(zeroIndexes) = [];
    % denominators(zeroIndexes) = [];
    
    % Now get the ratios of integers.
    ratios = numerators ./ (ones(s,1)*denominators);
    differences = abs(ratios - A_not_Good*ones(1,p));

    % Find the min difference:
    [minDifference, indexOfMin] = min(differences,[],2);

    partial_good = minDifference < args.tolerance;
    new_good = good;
    new_good(~good) = partial_good;
    
    numerators = diag(numerators(:,indexOfMin));
    denominators = denominators(:,indexOfMin);
    
    A_num(~good & new_good) = numerators(partial_good);
    A_den(~good & new_good) = denominators(partial_good);
    
    good = new_good;
       
    if prod(good) > 0
        break
    end
    
    min_denominator = max_denominator;
    max_denominator = max_denominator*args.step;

end

A_int = r(A_int);
A_num = r(A_num);
A_den = r(A_den);

A_sym = r(sym(A_num) ./ sym(A_den) + sym(A_int));
A_diff = abs(A_num ./ A_den + A_int - A_);
A_not_Good = prod(good > 0);

A_int = prepret;


    % 2019.11.20. (november 20, szerda), 14:32
    function s = prepret
        if args.structout
            s = struct();
            s.int = A_int;
            s.num = A_num;
            s.den = A_den;
            s.sym = A_sym;
            s.diff = A_diff;
            s.err = norm(A_diff);
            s.good = prod(good > 0);
            s.finalmaxden = max_denominator;
        else
            s = A_int;
        end
    end

end


function test1
%%

A = [
    1.00000000000001 1.34343434 3.411411411 12.111111111
    3.12312412313213 5.98732498 6.982199982 4.333333333
    ];

A_sym_builtin = sym(A)

A_sym = pcz_find_recdec(A)

end