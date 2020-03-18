function [ret] = pcz_feasible_PDLMI(PDLMI, plim, varargin)
%% pcz_feasible_PDLMI
%  
%  File: pcz_feasible_PDLMI.m
%  Directory: 2_demonstrations/lib/matlab
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  
%  Created on 2019. February 14.
%

%%

opts.title = 'PDLMI check';
opts.postol = 1e-6;
opts.tolerance = 1e-10;
opts.RandomPoints = 0;
opts.VectorArg = 1;
opts = parsepropval(opts, varargin{:});

TMP_PwVjnhhTEChSGsAudmgj = pcz_dispFunctionName(opts.title,'',struct('parent',1));

%%

PDLMI = value(PDLMI);

if ~any(isnan(PDLMI.Theta(:)) | isinf(PDLMI.Theta(:)))

    pcz_posdef_report_fh(PDLMI,PDLMI.subsvars,plim,opts);
        
else
    
    pcz_info(false, 'The solution contains NaN or Inf. The LMI is NOT feasible.')
    Coefficient_matrix = PDLMI.Theta;
    pcz_display(Coefficient_matrix);
    
end
    
pcz_dispFunctionEnd(TMP_PwVjnhhTEChSGsAudmgj);


end