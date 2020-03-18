function v = findsym(S,n)
%FINDSYM Finds the symbolic variables in a symbolic expression or matrix.
%
%   ================================================================
%   FINDSYM will be removed in a future release. Use SYMVAR instead.
%   ================================================================
%
%   FINDSYM(S), where S is a scalar or matrix sym, returns a string 
%   containing all of the symbolic variables appearing in S. The 
%   variables are returned in lexicographical order and are separated by
%   commas. If no symbolic variables are found, FINDSYM returns the
%   empty string.  The constants pi, i and j are not considered variables.
%
%   FINDSYM(S,N) returns the N symbolic variables closest to 'x' or 'X'. 
%   Upper-case variables are returned ahead of lower-case variables.
%   If S is a symbolic function the inputs to S are listed in front of the
%   other free variables.
%
%   Examples:
%      findsym(alpha+a+b) returns
%       a, alpha, b
%
%      findsym(cos(alpha)*b*x1 + 14*y,2) returns
%       x1, y
%
%      findsym(y*(4+3*i) + 6*j) returns
%       y

%   Copyright 1993-2014 The MathWorks, Inc.

Sym = privResolveArgs(S);
S = Sym{1};
if nargin == 2
    validateattributes(n,{'double'},{'scalar','positive','finite','integer'},'','N',2);
    scell = {};
    if isa(S,'symfun')
        argcell = privToCell(argnames(S));
        scell = cellfun(@(x){x.s}, argcell);
    end
    v = mupadmex('symobj::findsym', S.s, num2str(n), scell{:}, 0);
else
    v = mupadmex('symobj::findsym', S.s, 0);
end
v(v==' ')=[];
v = v(2:end-1); % trim quotes from around output
