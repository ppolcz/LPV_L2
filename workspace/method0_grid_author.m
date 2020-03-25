function method0_grid_author(modelname, A_fh, B_fh, C, D, p_lim, dp_lim, Resolution)
%%

pargrd = [ cellfun(@(c) {linspace(c{:})}, num2cell(num2cell([p_lim Resolution']),2))' num2cell(dp_lim,2)'];

% P_generate_symvars(1,p_lim(1),1,1);


p = sym('p',[3 1]);
dp = sym('dp',[3 1]);
pdp = [ p ; dp ];

p1 = p(1);
p2 = p(2);
p3 = p(3);

base = [
                   1
                  p1
                  p2
                  p3
 -(p1*p2)/(p2^2 - 7)
      -p1/(p2^2 - 7)
      -p2/(p2^2 - 7)
     p2^2/(p2^2 - 7)
 -(p2*p3)/(p2^2 - 7)
      -p3/(p2^2 - 7)
    ];
  
dbase = jacobian(base,p)*dp;

Pbase = matlabFunction(base, 'vars', {pdp.'});
Pbase_der = matlabFunction(dbase, 'vars', {pdp.'});

OVERALL = tic;
[gamma,Pvars] = lpvL2gain(@lpvsys,pargrd,Pbase,Pbase_der);

info = sprintf('(%s)x(%s)', ...
    strjoin(cellfun(@(n) {num2str(n)}, num2cell(Resolution) ), 'x'),...
    strjoin(cellfun(@(n) {num2str(n)}, num2cell(0*Resolution + 2) ), 'x'));

store_results('Results_All.csv', modelname, 0, gamma, toc(OVERALL), ...
    info, 'Tamas - grid')

function [Ap,Bp,Cp,Dp,S]=lpvsys(p)
    Ap = A_fh(p(1), p(2), p(3));
    Bp = B_fh(p(1), p(2), p(3));
    Cp = C;
    Dp = D;
end

end