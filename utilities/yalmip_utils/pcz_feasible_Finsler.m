function [ret] = pcz_feasible_Finsler(P, L, N, x_lim, varargin)
%% pcz_feasible_Finsler
%  
%  File: pcz_feasible_Finsler.m
%  Directory: 2_demonstrations/lib/matlab
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2018. May 22.
%

%%

opts.title = 'Finsler''s lemma check';
opts.postol = 1e-6;
opts.tolerance = 1e-10;
opts = parsepropval(opts, varargin{:});

TMP_PwVjnhhTEChSGsAudmgj = pcz_dispFunctionName(opts.title);

%%

Mode_xlim = 0;

PisAffinMatrix = isa(P,'PAffineMatrix');


if size(x_lim,2) == 2 && all(x_lim * [-1 ; 1] > 0) && numel(N.subsvars) == size(x_lim,1)
    [X_v,~,~] = P_ndnorms_of_X(x_lim);
    Mode_xlim = 1;
else
    X_v = x_lim;
end

if Mode_xlim
    try
        Mode_xlim = size(L,2) == size(N(X_v(1,:)'),1);
    catch
        Mode_xlim = 0;
    end
    
    if ~Mode_xlim
        X_v = x_lim;
    end
end


X_Nr = size(X_v,1);
P = value(P);
L = value(L);

Mode_str = { 'vertices are given' 'rectangular region: limits are given' };
pcz_info('Mode: %s, nr. of corners: %d', Mode_str{Mode_xlim+1}, X_Nr);
pcz_info('Tolerance: %g, positive tolerance: %g.', opts.tolerance, opts.postol)
pcz_info('LMI size: (%dx%d), annihilator rows: %d. ', size(P), size(L,2));
% pcz_dispFunction2('X_v = %s', pcz_num2str(X_v));
% pcz_dispFunction ''


feasible = zeros(X_Nr,1);
min_eig = zeros(X_Nr,1);
max_eig = zeros(X_Nr,1);
percentage = zeros(X_Nr,1);
iszero = zeros(X_Nr,1);


if PisAffinMatrix
    PL = [ P.Theta L ];
else
    PL = [ P L ];
end

PL_nan = isnan(PL);
PL_inf = isinf(PL);

Is_a_number = all(~PL_nan(:)) && all(~PL_inf(:));

if Is_a_number 

    for i = 1:X_Nr
        
        x_num = X_v(i,:)';
        
        if PisAffinMatrix
            P_num = P(x_num);
        else
            P_num = P;
        end
        
        N_num = N(x_num);
        
        % pcz_display(P_num,L,N_num,L*N_num,N_num'*L')

        [eigvals,min_eig(i),percentage(i)] = pcz_posdef_report(P_num + L*N_num + N_num'*L', opts.tolerance);
        feasible(i) = min_eig(i) > -opts.tolerance;
        max_eig(i) = max(eigvals);
        
        if max_eig(i) < opts.postol
            iszero(i) = 1;
        end
    end

    bool = all(feasible) && ~any(iszero);

    for i = 1:X_Nr
        if min_eig(i) < 0 && feasible(i)
            pcz_warning(false, 'Least eigv: %g, in corner %s', min_eig(i), pcz_num2str(X_v(i,:)));
        elseif ~feasible(i)
            pcz_info(false, 'Least eigv: %g, in corner %s', min_eig(i), pcz_num2str(X_v(i,:)));
        end
        
        if iszero(i)
            pcz_info(false, 'P+LN+N''L'' is almost zero (max eig: %d) in corner %s', max_eig(i), pcz_num2str(X_v(i,:)));
        end
    end

    pcz_dispFunctionSeparator
    
    if any(iszero)
        [indices,~,overallperc] = pcz_symzero_report(iszero);  
        pcz_info(false, 'P+LN+N''L'' is almost zero in %d corner points out of %d, (%g%%).', numel(indices), X_Nr, overallperc);
    end
    
    if ~bool
        [indices,~,overallperc] = pcz_symzero_report(~feasible);    
        pcz_info(false, 'This LMI is NOT feasible in %d corner points out of %d, (%g%%).', numel(indices), X_Nr, overallperc*100);
    else
        pcz_info(true, 'This LMI is feasible along the given tolerance value.');
    end
        
else
    
    pcz_info(false, 'Matrix P or L contains NaN of Inf. The LMI is NOT feasible.')
    
end
    
pcz_dispFunctionEnd(TMP_PwVjnhhTEChSGsAudmgj);


end
%%


