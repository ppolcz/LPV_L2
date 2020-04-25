function [ret] = P_generate_symvars_v10(nx, np, nw, nz, varargin)
%% 
%  
%  file:   P_generate_symvars_v5.m
%  author: Polcz Péter <ppolcz@gmail.com> 
%  
%  Created on 2017.02.13. Monpay, 14:49:31
%

%% beginning of the scope
% TMP_wtNPjzxHKNoJIigzXrEl = pcz_dispFunctionName;

%%

opts.verbose = 0;
opts = parsepropval(opts, varargin{:});

syms t real

if nargin <= 1, np = 0; end
if nargin <= 2, nw = 0; end
if nargin <= 3, nz = nw; end
if nargin <= 4, opts.verbose = 0; end

cmds = {
    't = sym(''t'',''real'');'
    pcz_generateSymStateVector(nx,'x','real')
    pcz_generateSymStateVector(np,'p','real')
    pcz_generateSymStateVector(np,'dp','real')
    pcz_generateSymStateVector(nw,'w','real')
    pcz_generateSymStateVector(nz,'z','real')
    
    'xw = [x ; w];'
    'xp = [x ; p];'
    'wp = [w ; p];'
    'xwp = [x ; w ; p];'
    'pdp = [p ; dp];'
    'xpdp = [x ; p ; dp];'
    'xwpdp = [x ; w ; p ; dp];'
    
    'nx = x_n;'
    'np = p_n;'
    'nw = w_n;'
    'nz = z_n;'
    'nxp = nx + np;'
    'nxp = nx + np;'
    'nwp = nw + np;'
    'nxpdp = nxp + np;'
    };

for i = 1:numel(cmds)
    cmd = cmds{i};
    if opts.verbose
        disp(cmd)
    end
    evalin('caller', cmd);
end

%% end of the scope
% pcz_dispFunctionend(TMP_wtNPjzxHKNoJIigzXrEl);

end