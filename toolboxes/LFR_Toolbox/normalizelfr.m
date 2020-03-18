% NORMALIZELFR - Normalisation of lfr-object
%-----------------------------------------------------------------
% PURPOSE
% Normalizes  the  variations  of REAL uncertain  parameters of an
% lfr-object  (it  also  applies  to  complex uncertain parameters
% with REAL nominal value).
%
% SYNOPSIS
% sys_out = normalizelfr(sys);
% sys_out = normalizelfr(sys,par_name);
% sys_out = normalizelfr(sys,par_name,d1,d2[,d3]);
%
% DESCRIPTION
% Let delta_i be the ith uncertain parameter.
% 1- min/max-type bounds: d1(i) < delta_i < d2(i). Calling
% 'normalizelfr' returns 'sys_out' such that
%       delta_i := (d2(i) + d1(i))/2 + delta_i * (d2(i) - d1(i))/2
% 2- min/max-type bounds with specified nomunal value. If the
% nominal value is not equal to (d2(i) + d1(i))/2:
%        delta_i := F(delta_i) where
%        F(-1) = d1(i), F(1) = d2(i), F(0) = d3(i)
% 3- disc-type bounds: |d1(i) - delta_i| < d2(i). Calling
% 'normalizelfr' returns 'sys_out' such that
%       delta_i := d1(i) + delta_i * d2(i)
%
% Note,  if d1, d2 and d3 are explicitely given, the corresponding
% values of the lfr-object are overwritten.
%
% NOTE
% If d1(i) = d2(i)  for  min/max-type  bounds  or if d2(i) = 0 for
% disc-type  bounds the uncertain parameter delta_i is substituted
% by  its  nominal value and vanishes from the uncertainty matrix.
%
% INPUT ARGUMENTS
% sys        lfr-object.
% par_name   Array  of  strings:  names  of  the  parameters to be
%            normalized.
% d1,d2      Vectors  of s ame length as par_name. The entries are
%            real  (bounds  and  radius)  or complex (disc center)
%            depending on the corresponding bound type.
% d3         Nominal value if not centered, only for real min/max-
%            type bounds (values ignored for disc-type bounds).
%
% OUTPUT ARGUMENT
% sys_out    Returned normalized lfr-object.
%
% See also lfr
%#----------------------------------------------------------------
% % EXAMPLE
% % Example 1: min/max-values of the lfr-object are overwritten
%    lfrs par1 par2 par3
%    sys1 = par1+1/par1+par2+par3;
%    par_name={'par1','par2','par3'};
%    dmin = [2; 0.1; 3];
%    dmax = [6; 0.2; 3];
%    sysn1 = normalizelfr(sys1,par_name,dmin,dmax)
% % Same result in a different way
%    lfrs par1 par2 par3 [2 0.1 3] [6 0.2 3]
%    sys2 = par1+1/par1+par2+par3;
%    sysn2 = normalizelfr(sys2);
%    distlfr(sysn1,sysn2)
%
% % Example 2: evaluation before and after normalization
%    sys = bnds2lfr('Y',[1 2;3 4],10,'p');
%    size(sys)
%    uplft(sys,{'Y1_1','Y1_2','Y2_1','Y2_2'},[0.9,1.8,3.3,3.6])
% % Look for same result after normalization
%    sysn = normalizelfr(sys);
%    size(sysn)
%    uplft(sysn,{'Y1_1','Y1_2','Y2_1','Y2_2'},[-1,-1,1,-1])
