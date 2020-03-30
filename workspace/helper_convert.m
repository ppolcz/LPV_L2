function [A_,B_,C_,D_,AC_,BD_,M_,nx,np,nu,ny] = helper_convert(A, B, C, D, arg5, arg6, arg7)
%% 
%  File: helper_convert.m
%  Directory: 8_published/LPV_L2/workspace
%  Author: Peter Polcz (ppolcz@gmail.com)
%
%  Created on 2020. March 25. (2019b)
%  Major review on 2020. March 30. (2019b)
% 
% Supported conversions:
% 
%  (A_fh,B_fh,C_fh,D_fh,p_lim), Class = 'lfr'
%  (A_fh,B_fh,C_fh,D_fh,p_lim,Class)
%  (A_fh,B_fh,C_fh,D_fh,p_lim,dp_lim), Class = 'lfr'
%  (A_fh,B_fh,C_fh,D_fh,p_lim,dp_lim,Class)
%  (A_lfr,B_lfr,C_lfr,D_lfr,Class)
% 
%  lfr|pfrl -> umat     (ureal
%              pmat     (pgrid)
%              plftmat  (tvreal)
% 
%  fh       -> lfr
%              umat
%              pmat
%              plftmat


%% Parse input arguments

p_lim = [];
dp_lim = [];

if nargin == 5 && ischar(arg5) % (A_lfr,B_lfr,C_lfr,D_lfr,Class)
    % A,B,C,D must be `lfr' or `plfr'
    
    A = plfr(A);
    B = plfr(B);
    C = plfr(C);
    D = plfr(D);
    
    Src_Class = 'plfr';
    Dst_Class = arg5;
    
elseif nargin == 5 % (A_fh,B_fh,C_fh,D_fh,p_lim), Class = 'lfr'
    % A,B,C,D must be `fh(p1,p2,...)' [default operation mode]

    Src_Class = 'fh';
    Dst_Class = 'lfr';

    p_lim = arg5;

elseif nargin == 6 && ischar(arg6) % (A_fh,B_fh,C_fh,D_fh,p_lim,Class)
    % A,B,C,D must be `fh(p1,p2,...)'

    Src_Class = 'fh';
    Dst_Class = arg6;

    p_lim = arg5;
    
elseif nargin == 6 % (A_fh,B_fh,C_fh,D_fh,p_lim,dp_lim), Class = 'lfr'
    % A,B,C,D must be `fh(p1,p2,...)' [default operation mode with pd_lim]
    
    Src_Class = 'fh';
    Dst_Class = 'lfr';
    
    p_lim = arg5;
    dp_lim = arg6;
  
elseif nargin == 7 % (A_fh,B_fh,C_fh,D_fh,p_lim,dp_lim,Class)
    % A,B,C,D must be `fh(p1,p2,...)'
    
    Src_Class = 'fh';
    Dst_Class = arg7;
    
    p_lim = arg5;
    dp_lim = arg6;

else
    error('Not supported')
end

if ~isa(p_lim,'cell')
    [~,p_cell] = pcz_generateLFRStateVector('p',p_lim,dp_lim);
else
    p_cell = p_lim;
end

%%

Cast = @(M) sprintf('%s(%s)',Dst_Class,inputname(1));

A_ = convert(A);
B_ = convert(B);
C_ = convert(C);
D_ = convert(D);

AC_ = [
    A_
    C_
    ];

BD_ = [
    B_
    D_
    ];

M_ = [ AC_ BD_ ];

np = numel(p_cell);
nx = size(A_,1);
nu = size(B_,2);
ny = size(C_,1);

    function dstobj = convert(srcobj)

        if isdouble(srcobj)
            
            % double -> CLASS
            dstobj = eval(Cast(srcobj));
            
        elseif strcmp(Dst_Class,'lfr') && ismember(class(srcobj),{'lfr','plfr'})
        
            % lfr|plfr -> lfr
            dstobj = lfr(srcobj);
            
        elseif ismember(Dst_Class,{'plftmat','umat'}) && ismember(class(srcobj),{'lfr','plfr'})
            
            % lfr|plfr -> plfr -> plftmat|ureal
            srcobj = plfr(srcobj);
            dstobj = eval(Cast(srcobj));
            
        else
            
            % lfr -> plfr
            if isa(srcobj,'lfr')
                srcobj = plfr(srcobj);
            end
            
            % fh|plfr -> CLASS
            dstobj = srcobj(p_cell{:});
        end
    end

end