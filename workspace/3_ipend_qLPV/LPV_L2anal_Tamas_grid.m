function LPV_L2anal_Tamas_grid(A_fh, B_fh, C, D, K, x_lim, p_lim, dp_lim, LPVTools, varargin)
%%

pargrd = [ cellfun(@(c) {linspace(c{:})}, num2cell(num2cell([p_lim LPVTools.Resolution]),2))' num2cell(dp_lim,2)'];

Pbase = @(p) [
    1 
    p(1) 
    p(2) 
    p(3) 
%     1/(p(2)^2 - 2)
    ];
Pbase_der = @(p) [
    0 
    p(4) 
    p(5) 
    p(6)                
%     -2*p(2)*p(5) / (p(2)^2 - 2)^2
    ];     

[gam,Pvars] = lpvL2gain(@lpvsys,pargrd,Pbase,Pbase_der) 

function [Ap,Bp,Cp,Dp,S]=lpvsys(p)
    Ap = A_fh(p(1), p(2), p(3));
    Bp = B_fh(p(1), p(2), p(3));
    Cp = C;
    Dp = D;
end

end