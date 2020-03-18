function [ret] = pcz_symrac(val,dim)
%% Script pcz_symrac
%  
%  File: pcz_symrac.m
%  Directory: /1_projects/2_sta/powertrain
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2018. March 17.
%

%%

if nargin < 2
    dim = 2;
end

if isscalar(val)

    ret = struct;
    
    [num,den] = numden(val);
    
    ret.ispoly = isempty(symvar(den));
        
    [ret.b,ret.tb] = coeffs(num);
    [ret.a,ret.ta] = coeffs(den);
    
    ret.b = ret.b / ret.a(1);
    ret.a = ret.a / ret.a(1);

    if ret.ispoly && ret.a ~= 1
        ret.a
        error 'JAVITSD KI!!! ENNEK EGYNEK KELL LENNIE!'
    end

    ret.b = double(ret.b);
    ret.a = double(ret.a);
    
    return
end


ret = 'TODO';

end