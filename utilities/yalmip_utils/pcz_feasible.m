function [ret] = pcz_feasible(sol, CONS, varargin)
%% pcz_feasible
%  
%  File: pcz_feasible.m
%  Directory: 2_demonstrations/lib/matlab
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2018. May 18.
%

%%

opts.tolerance = 1e-10;
opts = parsepropval(opts, varargin{:});

%%

pcz_info(sol.problem == 0, '%s. Solver time: <strong>%g</strong>', ...
    sol.info, sol.solvertime);

[Prim,Dual] = check(CONS);

mp = min(Prim);
md = min(Dual);

msg_tolerance = sprintf('Tolerance: %g.', opts.tolerance);

msg1 = [ 'The solution is feasible. ' msg_tolerance ];
msg2 = [ 'The solution is NOT feasible. ' msg_tolerance ];

msg3 = '';
if mp <= 0
    msg3 = sprintf(' Min(Primal) = %d.', mp);
end

if md <= 0
    msg3 = sprintf(' Min(Dual) = %d.', md);
end

bool = mp < -opts.tolerance || md < -opts.tolerance;

if bool
    pcz_info(false, [msg2 msg3]);
else
    pcz_warning(mp > 0 && md > 0, [msg1 msg3]);
end

if nargout > 0
    ret = bool;
end

end