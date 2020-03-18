% UPLFT        - Upper linear fractional transformation
%-----------------------------------------------------------------
% PURPOSE
% Partial  or  full  upper  linear  fractional transformation. 1/s
% equal to infinity is supported.
%
% SYNOPSIS
% [sysout,iserror] = uplft(sysin,names,values);
% [sysout,iserror] = uplft(sysin,names,values,names,values,...);
%
% INPUT ARGUMENTS
% sysin     lfr-object
% names     cells  or character arrays e.g., {'Int','b'}, {Int,b},
%           {'1/s','b'} or char('Int','b').
%           The  entries  of  the cell can be strings, symbolic or
%           lfr elementary objects.
%           In  case  of duplicate names, it is the last occurence
%           that is considered.
% values    row-vector or row-cell.
%           Entries with values equal to NaN are ignored.
%           For  DC-gain  computation, set values{ii} = Inf for ii
%           s.t. names{ii} = 'Int'.
%           If  values{ii} is a scalar and the corresponding block
%           is  full, this block is replaced by a random blocks of
%           magnitude given by the scalar value
%
% OUTPUT ARGUMENTS
% sysout    sysin after partial Delta-loop is closed
% iserror   iserror = 1  if the result is not well-posed.
%
% See also lfr/eval, usubs
%#----------------------------------------------------------------
% % EXAMPLE
%    lfrs a b c
%    lfrobj = [1/a b^2 a+b+c;1+a^2 c*2*j/(a+b) 4];
%    x1 = uplft(lfrobj,{a,b,c},{1,2,3})
%    x2 = uplft(lfrobj,'a',1,'b',2,'c',3)
%    x3 = uplft(lfrobj,char('a','b','c'),{1,2,3})
%
% % Partial upper LFT
%    x4 = uplft(lfrobj,{a,c},{1,3})
%
% % DC-gain
%    lfrs Int
%    lfrobj = [1/(1+a*Int) b^2 a+b+c;1/(1+a^2*Int) c*2*j/(1+a+b) 4];
%    x5 = uplft(lfrobj,{Int,a,b,c},{Inf,3,2,3})
%    x6 = uplft(lfrobj,{Int,a,b,c},{1e+8,3,2,3})
%    x5.d - x6.d
%
% % Full blocks
%    d = lfr('d','ltifc',[2,3]);
%    lfrobj = [ [a;b] d/b];
%    x7 = uplft(lfrobj,{a,b,d},{1,2,ones(2,3)})
%    x8 = uplft(lfrobj,{a,b,d},{1,2,3})
